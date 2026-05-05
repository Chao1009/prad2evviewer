#!/usr/bin/env python3
"""
post_elog.py — manually post a saved local auto-report to the JLab logbook.

Mirrors the prad2_server logic:

  * Dedup check: GET <elog_url>/api/elog/entries?book=<book>&title=...&
      field=lognumber&field=title&limit=20      (server: handleElogCheck)
    Substring server-side, EXACT title match client-side.

  * Post:        PUT <elog_url>/incoming/prad2_report.xml with the XML
                 body, cert + key auth                 (server: handleElogPost)

Stdlib only — no `requests`, no `pip install`.

Usage:
  ./post_elog.py 24364                       # latest XML in run_024364/
  ./post_elog.py 024364                      # leading zeros OK
  ./post_elog.py path/to/report.xml          # explicit path
  ./post_elog.py 24364 --check-only          # only print dedup result
  ./post_elog.py 24364 --force               # skip dedup check
"""

import argparse
import json
import re
import ssl
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

# Defaults match clonfarm11 / clonpc19 setup. Override via flags or by
# pointing --config at a monitor_config.json.
DEFAULT_REPORTS_DIR = "/home/clasrun/prad2_daq/monitor/reports"
DEFAULT_CERT        = "/home/clasrun/prad2_daq/monitor/keys/elog-cert.pem"
DEFAULT_KEY         = "/home/clasrun/prad2_daq/monitor/keys/elog-key.pem"
DEFAULT_URL         = "https://logbooks.jlab.org"
DEFAULT_BOOK        = "PRADLOG"
TITLE_BASE          = "PRad2 Event Monitor Auto Report"
TIMEOUT_S           = 30

# Tag substitutions applied to every upload so old saved XMLs (whose
# <Tags> block still references tags from before the elog enum was
# tightened) post cleanly.  No-op on new XMLs that already use the
# correct tags. Override / extend with --rewrite-tag OLD=NEW; disable
# entirely with --no-rewrite-tags.
DEFAULT_TAG_REWRITES = {
    "AutoReport": "Autolog",
    "PRad2":      "DAQ",
}


def auto_title(run: int) -> str:
    return f"Run #{run}: {TITLE_BASE}"


def find_xml_for_run(reports_dir: str, run: int):
    rundir = Path(reports_dir) / f"run_{run:06d}"
    if not rundir.is_dir():
        return None
    xmls = sorted(rundir.glob("report_*.xml"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return xmls[0] if xmls else None


def run_from_path(path):
    """Run number embedded in a path (e.g. .../run_024364/...), or None."""
    m = re.search(r"run_0*(\d+)", str(path))
    return int(m.group(1)) if m else None


def make_ssl_ctx(cert: str, key: str) -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.load_cert_chain(certfile=cert, keyfile=key)
    return ctx


def check_elog(url: str, book: str, cert: str, key: str, title: str):
    """Mirror handleElogCheck — exact-title match against the live logbook."""
    qs = urllib.parse.urlencode([
        ("book",  book),
        ("title", title),
        ("field", "lognumber"),
        ("field", "title"),
        ("limit", "20"),
    ])
    full = f"{url}/api/elog/entries?{qs}"
    ctx = make_ssl_ctx(cert, key)
    req = urllib.request.Request(full, method="GET")
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=TIMEOUT_S) as r:
            body = r.read().decode("utf-8", errors="replace")
            code = r.status
    except urllib.error.HTTPError as e:
        return {"checked": False, "detail": f"HTTP {e.code}"}
    except Exception as e:
        return {"checked": False, "detail": f"network: {e}"}
    if code != 200:
        return {"checked": False, "detail": f"HTTP {code}"}
    try:
        j = json.loads(body)
    except json.JSONDecodeError as e:
        return {"checked": False, "detail": f"json: {e}"}
    if j.get("stat") != "ok":
        return {"checked": False,
                "detail": j.get("message", "non-ok stat")}
    entries = j.get("data", {}).get("entries", [])
    for e in entries:
        if e.get("title") == title:
            return {"checked": True, "exists": True,
                    "lognumber": e.get("lognumber"),
                    "matched_count": len(entries)}
    return {"checked": True, "exists": False,
            "matched_count": len(entries)}


def rewrite_tags(xml_text: str, mapping: dict) -> str:
    """Apply --rewrite-tag old=new replacements inside <Tags><tag>…
    blocks. mapping is {"AutoReport": "Autolog", ...}.  Match is
    exact-text inside a <tag>…</tag> element so we don't accidentally
    munge body content.
    """
    if not mapping:
        return xml_text
    def _sub(m):
        old = m.group(1)
        new = mapping.get(old, old)
        return f"<tag>{new}</tag>"
    return re.sub(r"<tag>([^<]*)</tag>", _sub, xml_text)


def post_elog(url: str, cert: str, key: str, xml_path, run: int,
              rewrite_map: dict = None):
    """Mirror handleElogPost — shells out to curl --upload-file, the
    same command prad2_server runs.  Apache on logbooks.jlab.org is
    picky about the exact request shape (Expect: 100-continue, header
    set, etc.); using curl directly avoids tripping over those
    differences from urllib's PUT.

    Filename embeds the run number so the elog's per-/incoming/-name
    one-shot policy doubles as a per-run dedup guard: only the first
    post for a given run lands; subsequent attempts (other monitors,
    replay scripts) hit "already processed" and are rejected before
    creating duplicates.

    --rewrite-tags lets us replay older saved XMLs whose <Tags> block
    still references tags that aren't in the elog's current
    enumeration (e.g. AutoReport → Autolog) without mutating the
    on-disk archive.
    """
    if run > 0:
        upload_name = f"prad2_run_{run:06d}.xml"
    else:
        upload_name = "prad2_" + Path(xml_path).name
    full = f"{url}/incoming/{upload_name}"
    marker = "___HTTP_CODE___"

    if rewrite_map:
        # Stage to a temp file with the substituted tags so curl can
        # still --upload-file a real path.
        import tempfile
        text = Path(xml_path).read_text(encoding="utf-8")
        text = rewrite_tags(text, rewrite_map)
        tf = tempfile.NamedTemporaryFile(mode="w", suffix=".xml",
                                         delete=False, encoding="utf-8")
        tf.write(text)
        tf.close()
        upload_src = tf.name
    else:
        upload_src = str(xml_path)

    cmd = [
        "curl", "-sS",
        "--cert", cert, "--key", key,
        "--upload-file", upload_src,
        full,
        "-w", f"\n{marker}%{{http_code}}",
        "-m", str(TIMEOUT_S),
    ]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=TIMEOUT_S + 5)
    finally:
        if rewrite_map:
            try: Path(upload_src).unlink()
            except Exception: pass
    if isinstance(p, Exception):
        return 0, str(p)
    out = (p.stdout or "") + (p.stderr or "")
    m = re.search(re.escape(marker) + r"(\d+)\s*$", out)
    if not m:
        return 0, out.strip() or "no HTTP code from curl"
    code = int(m.group(1))
    body = out[:m.start()].rstrip()
    return code, body


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target",
        help="run number (e.g. 24364 or 024364) or path to a report.xml")
    ap.add_argument("--reports-dir", default=DEFAULT_REPORTS_DIR)
    ap.add_argument("--cert", default=DEFAULT_CERT)
    ap.add_argument("--key",  default=DEFAULT_KEY)
    ap.add_argument("--url",  default=DEFAULT_URL)
    ap.add_argument("--book", default=DEFAULT_BOOK)
    ap.add_argument("--run",  type=int,
        help="override run number when target is an XML path")
    ap.add_argument("--force", action="store_true",
        help="skip the dedup check (post even if already in elog)")
    ap.add_argument("--check-only", action="store_true",
        help="print dedup result and exit; do not post")
    ap.add_argument("--rewrite-tag", action="append", default=[],
        metavar="OLD=NEW",
        help="extra tag rewrite to apply before upload (repeatable). "
             f"defaults already cover {DEFAULT_TAG_REWRITES}")
    ap.add_argument("--no-rewrite-tags", action="store_true",
        help="disable the default tag rewrites (post the XML verbatim)")
    args = ap.parse_args(argv)

    # Resolve XML path + run number
    if Path(args.target).is_file():
        xml_path = Path(args.target)
        run = args.run if args.run is not None else run_from_path(xml_path)
    else:
        try:
            run = args.run if args.run is not None else int(args.target.lstrip("0") or "0")
        except ValueError:
            sys.exit(f"could not parse run number from '{args.target}'")
        xml_path = find_xml_for_run(args.reports_dir, run)
        if xml_path is None:
            sys.exit(f"no XML found in {args.reports_dir}/run_{run:06d}")
    if not run:
        sys.exit("could not determine run number")

    title = auto_title(run)
    print(f"title : {title}")
    print(f"xml   : {xml_path}")

    if not args.force:
        print("checking elog...")
        c = check_elog(args.url, args.book, args.cert, args.key, title)
        if not c.get("checked"):
            print(f"  dedup unavailable: {c.get('detail')}")
            if not args.check_only:
                # Match the server's fail-closed policy on uncertain dedup.
                sys.exit("dedup couldn't run — refusing to post (re-run with --force to override)")
        elif c.get("exists"):
            print(f"  ALREADY IN ELOG: lognumber={c.get('lognumber')}  (skipping)")
            sys.exit(0)
        else:
            print(f"  not in elog ({c.get('matched_count')} fuzzy hits)")
        if args.check_only:
            return

    print("posting...")
    rewrite_map = {} if args.no_rewrite_tags else dict(DEFAULT_TAG_REWRITES)
    for spec in args.rewrite_tag:
        if "=" not in spec:
            sys.exit(f"--rewrite-tag expects OLD=NEW, got '{spec}'")
        old, new = spec.split("=", 1)
        rewrite_map[old.strip()] = new.strip()
    if rewrite_map:
        print(f"tag rewrites: {rewrite_map}")
    code, body = post_elog(args.url, args.cert, args.key, xml_path, run,
                           rewrite_map=rewrite_map or None)
    print(f"HTTP {code}")
    if body:
        print(body[:1000])
    sys.exit(0 if code in (200, 201) else 1)


if __name__ == "__main__":
    main()
