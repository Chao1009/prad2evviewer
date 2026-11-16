#!/usr/bin/env python3
"""
lms_pass.py - LMS amplitude pass, epoch reducers and HV/LMS segments (plan P14.3: D4, D8).

Called through the snake_scan.py dispatcher (convention add_parser/run):
  snake_scan.py lms FILE...                          per-file extraction -> output/snake_scan/lms_npz/lms_<run>_<split>.npz
        [--map PATH] [--cache-dir DIR] [-j N] [--force]       (map: derived map for runs <= 23849, DB map above)
  snake_scan.py lms --reduce segments                -> out/segments_23701.tsv          (D8, 9 local 23701 files)
  snake_scan.py lms --reduce scan [--npz23744 P]     -> out/lms_scan.tsv (23701 alone, primary), out/lms_scan_merged.tsv (D4 merge)
  snake_scan.py lms --reduce prod [--rt RT ...]      -> out/lms_prod_<rt>.tsv           (D4, local production splits)

lms_pass(evio, map, cache=None) -> npz path.  Keys (PbWO4 modules by name, `mod_names`, order of the pipeline):
  lms_I[n, nW]    single-peak WaveAnalyzer integral of every PbWO4 channel in the bit-24 events (nan: absent,
                  npeaks != 1 or integral <= 0)
  lms_ref[n, 3]   same for LMS1..3 (looked up by NAME through the map, never by slot)
  alpha_ref[m, 3] LMS1..3 integrals in bit-25 events without bit 24 (each alpha event fires about one PMT)
  lms_ts, alpha_ts  TI timestamps (ticks of 4 ns);   lms_nch[n] channels with samples (all detectors);
  lms_pk_med[n]   median over PbWO4 modules of the peak height above pedestal (ADC);
  ped[nW]         median over the LMS events of the WaveAnalyzer pedestal mean, per module (D8);  n_all  events
Normalisation (GainCorrCompute convention, as lms_gain_variation.py): per event
  I_norm,i = I_i x mean_{j=LMS2,LMS3}( A_j / ref_j ),   A_j = median alpha integral of the same file;
a file with fewer than MIN_ALPHA alpha events takes A_j from all alpha events of the reduced group (epoch).
L_i = median of I_norm,i over the group's LMS events, err = 1.25 x 1.48 MAD / sqrt(n).
"""
from __future__ import annotations

import ast
import glob
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
PYSCRIPTS = SS.parents[1]
OUT = SS / "out"
NPZ_DIR = PYSCRIPTS / "output" / "snake_scan" / "lms_npz"     # per-split cache, in the git-ignored output area
EVIO_DIR = Path("/Users/L00338853/work/PRad/data/evio")
BENCH = PYSCRIPTS / "waveform_analysis" / "energy_resolution_benchmark.py"
TBIT_LMS, TBIT_ALPHA = 1 << 24, 1 << 25
REF_USE = (1, 2)            # LMS2, LMS3 (0-based)
MIN_ALPHA = 10              # alpha events per file below which the group's alpha median is used
MIN_N = 3                   # events per module for a reported L
ERR_K = 1.25 * 1.4826
DERIVED_MAP_MAX_RUN = 23849
SEG_SIGMA = 3.0
for p in (str(SS), str(PYSCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)


def _run_split(path):
    m = re.search(r"prad_0*(\d+)\.evio\.0*(\d+)$", str(path))
    return int(m.group(1)), int(m.group(2))


def default_map(path):
    import snake_scan as SSN
    return str(SSN.DERIVED_MAP) if _run_split(path)[0] <= DERIVED_MAP_MAX_RUN else ""


def cache_path(path, cache_dir=NPZ_DIR):
    r, s = _run_split(path)
    return Path(cache_dir) / f"lms_{r}_{s:05d}.npz"


# ---------------------------------------------------------------------------------------------------------------
# extraction
# ---------------------------------------------------------------------------------------------------------------

def lms_pass(evio, map="", cache=None, force=False):
    """One EVIO file -> npz (see module docstring).  map "" = the database hycal_map.  Returns the npz path."""
    import _common as C
    from prad2py import dec
    evio = str(evio)
    cache = Path(cache) if cache else cache_path(evio)
    if cache.exists() and not force:
        return cache
    t0 = time.time()
    pipe = C.setup_pipeline(evio_path=evio, hc_map_file=str(map) if map else "")
    hy = pipe.hycal
    mods = [hy.module(i) for i in range(hy.module_count())]
    pw = [i for i, m in enumerate(mods) if m.is_pwo4()]
    col = {mi: k for k, mi in enumerate(pw)}
    ref_idx = {i: int(m.name[3]) - 1 for i, m in enumerate(mods)
               if len(m.name) == 4 and m.name.startswith("LMS") and m.name[3] in "123"}
    nW = len(pw)
    ch = dec.EvChannel()
    ch.set_config(pipe.cfg)
    if ch.open_auto(evio) != dec.Status.success:
        raise RuntimeError(f"cannot open {evio}")
    lms_ts, lms_I, lms_ref, lms_nch, lms_pk, lms_ped = [], [], [], [], [], []
    alp_ts, alp_ref = [], []
    n_all = 0
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            tb = int(inf.trigger_bits)
            n_all += 1
            is_l, is_a = bool(tb & TBIT_LMS), bool(tb & TBIT_ALPHA)
            if not (is_l or is_a):
                continue
            e = ch.fadc()
            row = np.full(nW, np.nan, np.float32)
            pk = np.full(nW, np.nan, np.float32)
            pd = np.full(nW, np.nan, np.float32)
            ref = np.full(3, np.nan, np.float32)
            nch = 0
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present:
                    continue
                crate = pipe.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        cd = slot.channel(c)
                        if cd.nsamples <= 0:
                            continue
                        nch += 1
                        m = hy.module_by_daq(crate, s, c)
                        if m is None:
                            continue
                        mi = m.index
                        if mi not in ref_idx and not (is_l and mi in col):
                            continue
                        pm, _, pks = pipe.wave_ana.analyze(cd.samples)
                        if mi in ref_idx:
                            if len(pks) == 1 and pks[0].integral > 0:
                                ref[ref_idx[mi]] = pks[0].integral
                        else:
                            pd[col[mi]] = pm
                            if len(pks) == 1 and pks[0].integral > 0:
                                row[col[mi]] = pks[0].integral
                                pk[col[mi]] = pks[0].height
            ts = int(inf.timestamp)
            if is_l:
                lms_ts.append(ts); lms_I.append(row); lms_ref.append(ref); lms_nch.append(nch)
                lms_pk.append(np.nanmedian(pk) if np.isfinite(pk).any() else np.nan); lms_ped.append(pd)
            else:
                alp_ts.append(ts); alp_ref.append(ref)
    ped = (np.nanmedian(np.asarray(lms_ped, np.float32), axis=0) if lms_ped
           else np.full(nW, np.nan, np.float32))
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache, mod_names=np.array([mods[i].name for i in pw]), n_all=n_all, source=os.path.basename(evio),
        map=os.path.basename(str(map)) if map else "hycal_map.json (database)",
        lms_I=np.asarray(lms_I, np.float32).reshape(-1, nW), lms_ref=np.asarray(lms_ref, np.float32).reshape(-1, 3),
        lms_ts=np.asarray(lms_ts, np.int64), lms_nch=np.asarray(lms_nch, np.int32),
        lms_pk_med=np.asarray(lms_pk, np.float32), alpha_ref=np.asarray(alp_ref, np.float32).reshape(-1, 3),
        alpha_ts=np.asarray(alp_ts, np.int64), ped=ped)
    print(f"[lms] {os.path.basename(evio)}: {n_all} events, {len(lms_ts)} LMS, {len(alp_ts)} alpha "
          f"({time.time() - t0:.0f} s)", flush=True)
    return cache


def ref_names_by_daq(map="", evio_path=""):
    """{(crate, slot, channel): module name} for LMS1-3 through the given map (test of the by-name rule)."""
    import _common as C
    pipe = C.setup_pipeline(evio_path=evio_path, hc_map_file=str(map) if map else "")
    hy = pipe.hycal
    out = {}
    for i in range(hy.module_count()):
        m = hy.module(i)
        if m.name in ("LMS1", "LMS2", "LMS3"):
            out[(int(m.daq.crate), int(m.daq.slot), int(m.daq.channel))] = m.name
    return out


# ---------------------------------------------------------------------------------------------------------------
# reducers
# ---------------------------------------------------------------------------------------------------------------

def load(p):
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files}


def _alpha_median(files):
    a = np.concatenate([f["alpha_ref"] for f in files if len(f["alpha_ref"])] or [np.empty((0, 3), np.float32)])
    return np.nanmedian(a, axis=0) if len(a) else np.full(3, np.nan)


def _n_alpha(f):
    """Alpha events with a finite integral, minimum over LMS2 and LMS3 (each alpha event fires about one PMT)."""
    return int(np.isfinite(f["alpha_ref"][:, list(REF_USE)]).sum(axis=0).min()) if len(f["alpha_ref"]) else 0


def normalise(f, a_group):
    """(I_norm[n, nW], used_group_alpha) for one file dict; A_j = file alpha median, else the group's."""
    used_group = _n_alpha(f) < MIN_ALPHA
    A = a_group if used_group else np.nanmedian(f["alpha_ref"], axis=0)
    r = f["lms_ref"][:, list(REF_USE)]
    fac = np.mean(A[list(REF_USE)][None, :] / r, axis=1)           # nan if either reference is missing
    return f["lms_I"] * fac[:, None], used_group


def median_err(x):
    """Per column median, 1.25 x 1.48 MAD / sqrt(n), n (nan ignored)."""
    n = np.isfinite(x).sum(axis=0)
    with np.errstate(all="ignore"):
        med = np.nanmedian(x, axis=0)
        mad = np.nanmedian(np.abs(x - med), axis=0)
    err = ERR_K * mad / np.sqrt(np.maximum(n, 1))
    ok = n >= MIN_N
    return np.where(ok, med, np.nan), np.where(ok, err, np.nan), n


def reduce_group(files):
    """files: list of loaded npz dicts of one epoch -> (names, L, err, n, info)."""
    names = [str(s) for s in files[0]["mod_names"]]
    a_group = _alpha_median(files)
    blocks, n_group = [], 0
    for f in files:
        if [str(s) for s in f["mod_names"]] != names:
            order = {str(s): k for k, s in enumerate(f["mod_names"])}
            idx = [order[nm] for nm in names]
            f = dict(f, lms_I=f["lms_I"][:, idx])
        x, ug = normalise(f, a_group)
        n_group += int(ug)
        blocks.append(x)
    L, err, n = median_err(np.concatenate(blocks, axis=0))
    info = dict(n_lms=int(sum(len(f["lms_ts"]) for f in files)), n_files=len(files), n_files_group_alpha=n_group,
                n_alpha=int(sum(len(f["alpha_ref"]) for f in files)))
    return names, L, err, n, info


def write_tsv(path, names, L, err, n, source, header_lines, extra_col=None, extra_name="source"):
    with open(path, "w") as fh:
        for h in header_lines:
            fh.write(f"# {h}\n")
        fh.write("# source: " + source + "\n")
        fh.write("module\tL\terr\tn" + (f"\t{extra_name}" if extra_col is not None else "") + "\n")
        for k, nm in enumerate(names):
            if not np.isfinite(L[k]):
                continue
            fh.write(f"{nm}\t{L[k]:.4f}\t{err[k]:.4f}\t{int(n[k])}")
            if extra_col is not None:
                fh.write(f"\t{extra_col[k]}")
            fh.write("\n")


def run_files(run):
    return sorted(NPZ_DIR.glob(f"lms_{run}_*.npz"))


def reduce_segments(run=23701, out=None):
    """D8: per-file LMS2/LMS3 alpha-normalised reference R_j = median(LMS_j) / median(alpha_j), median PbWO4
    pedestal; a step is a consecutive difference > 3 sigma (1.48 MAD of the differences of that series)."""
    paths = run_files(run)
    files = [load(p) for p in paths]
    A = _alpha_median(files)
    rows, series = [], {"R2": [], "R3": [], "ped": []}
    for p, f in zip(paths, files):
        Aj = np.nanmedian(f["alpha_ref"], axis=0) if _n_alpha(f) >= MIN_ALPHA else A
        R = np.nanmedian(f["lms_ref"], axis=0) / Aj
        pm = float(np.nanmedian(f["ped"]))
        series["R2"].append(R[1]); series["R3"].append(R[2]); series["ped"].append(pm)
        rows.append((p.stem.split("_")[-1], len(f["lms_ts"]), len(f["alpha_ref"]), int(f["lms_ts"][0]) if len(f["lms_ts"]) else 0,
                     R[1], R[2], pm))
    steps, zs = [], {}
    for k, v in series.items():
        v = np.asarray(v, float)
        ok = np.where(np.isfinite(v))[0]          # files without LMS events (nan) are skipped in the chain
        d = np.diff(v[ok])
        sig = ERR_K / 1.25 * np.nanmedian(np.abs(d - np.nanmedian(d)))
        z = d / sig if sig > 0 else np.full(len(d), np.nan)
        zs[k] = (ok, d, z, sig)
        for j in np.where(np.abs(z) > SEG_SIGMA)[0]:
            steps.append((k, rows[ok[j]][0], rows[ok[j + 1]][0], d[j], z[j]))
    out = Path(out) if out else OUT / f"segments_{run}.tsv"
    with open(out, "w") as fh:
        fh.write("# D8 segments, run %d: per-file LMS2/LMS3 (median LMS / median alpha) and median PbWO4 pedestal (ADC)\n" % run)
        fh.write("# source: " + " ".join(p.name for p in paths) + "; map hycal_map_til_run023849.json\n")
        fh.write("# step = consecutive-file difference > %.0f sigma, sigma = 1.48 MAD of that series' differences\n" % SEG_SIGMA)
        for k in series:
            fh.write(f"# sigma_{k} = {zs[k][3]:.5g}\n")
        fh.write("split\tn_lms\tn_alpha\tts_first\tR_LMS2\tR_LMS3\tped_med\tz_R2\tz_R3\tz_ped\n")
        for j, r in enumerate(rows):
            zz = []
            for k in ("R2", "R3", "ped"):          # z of the step from the previous valid file to this one
                ok = list(zs[k][0])
                zz.append("%.2f" % zs[k][2][ok.index(j) - 1] if j in ok and ok.index(j) > 0 else "-")
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.5g}\t{r[5]:.5g}\t{r[6]:.5g}\t" + "\t".join(zz) + "\n")
        if steps:
            fh.write("# steps: " + "; ".join(f"{k} {a}->{b} d={d:.4g} z={z:.1f}" for k, a, b, d, z in steps) + "\n")
        else:
            fh.write("# steps: none (no |z| > %.0f)\n" % SEG_SIGMA)
    return steps, rows, zs


def reduce_scan(npz_23744=None, out=None):
    """D4 epoch values for the scan.  Primary (out, default lms_scan.tsv): 23701 alone (D4 "else 23701 alone"; 23744
    sits about 2.5 % below 23701, so a per-module merge would turn that common factor into scatter).  Variant
    (<out stem>_merged.tsv, written when the 23744 npz exists): 23701 merged with 23744 where it agrees within 2 SE.
    Returns the primary arrays; stats carries the merged ones (L_merged, err_merged, n_merged, src_merged)."""
    paths = run_files(23701)
    f01 = [load(p) for p in paths]
    names, L1, e1, n1, i1 = reduce_group(f01)
    L, err, n = L1.copy(), e1.copy(), n1.copy()
    src = np.array(["23701"] * len(names), dtype=object)
    Lp, errp, np_, srcp = L1.copy(), e1.copy(), n1.copy(), src.copy()      # primary: 23701 alone
    hdr = [f"23701: {i1['n_lms']} LMS events in {i1['n_files']} files, {i1['n_alpha']} alpha events; "
           f"{i1['n_files_group_alpha']} files used the epoch alpha median (< {MIN_ALPHA} alpha events)"]
    srcline = " ".join(p.name for p in paths) + "; map hycal_map_til_run023849.json"
    stats = dict(i1=i1)
    out = Path(out) if out else OUT / "lms_scan.tsv"
    merged_out = out.with_name(out.stem + "_merged.tsv")
    if npz_23744 and Path(npz_23744).exists():
        f44 = load(npz_23744)
        if "alpha_ref" not in f44 or not len(f44["alpha_ref"]):
            f44["alpha_ref"] = np.empty((0, 3), np.float32)
        if "lms_ts" not in f44:
            f44["lms_ts"] = f44.get("ts", np.zeros(len(f44["lms_I"]), np.int64))
        names2, L2, e2, n2, i2 = reduce_group([f44]) if _n_alpha(f44) >= MIN_ALPHA else _reduce_with(f44, f01, names)
        idx = {nm: k for k, nm in enumerate(names2)}
        j = np.array([idx.get(nm, -1) for nm in names])
        ok = (j >= 0)
        L2m = np.where(ok, L2[np.maximum(j, 0)], np.nan); e2m = np.where(ok, e2[np.maximum(j, 0)], np.nan)
        n2m = np.where(ok, n2[np.maximum(j, 0)], 0)
        both = np.isfinite(L1) & np.isfinite(L2m)
        with np.errstate(all="ignore"):
            z = (L1 - L2m) / np.hypot(e1, e2m)
        agree = both & (np.abs(z) <= 2.0)
        dis = both & ~agree
        w1, w2 = 1 / e1 ** 2, 1 / e2m ** 2
        with np.errstate(all="ignore"):
            Lm = (w1 * L1 + w2 * L2m) / (w1 + w2)
            em = 1 / np.sqrt(w1 + w2)
        L = np.where(agree, Lm, L); err = np.where(agree, em, err); n = np.where(agree, n1 + n2m, n)
        src = np.where(agree, "merged", src)
        stats.update(L_merged=L, err_merged=err, n_merged=n, src_merged=src, i2=i2, n_agree=int(agree.sum()),
                     n_dis=int(dis.sum()), n_both=int(both.sum()),
                     dis=[(names[k], float(L1[k]), float(L2m[k]), float(z[k])) for k in np.where(dis)[0]])
        hdr.append(f"23744: {i2['n_lms']} LMS events; merged {stats['n_agree']} of {stats['n_both']} modules "
                   f"(agree within 2 SE), 23701 alone for the rest; disagreements: {stats['n_dis']}")
        hdr.append("23744 sits about 2.5 % below 23701 (systematic, see disagreements); the merge turns that common factor into "
                   "scatter, so lms_scan.tsv (23701 alone) is the primary L_scan and this merged table is the sensitivity variant")
        for nm, a, b, zz in stats["dis"]:
            hdr.append(f"disagree {nm}: L23701={a:.4g} L23744={b:.4g} z={zz:.2f}")
        write_tsv(merged_out, names, L, err, n, srcline + " + " + os.path.basename(str(npz_23744)), hdr, src, "source")
        hdr = hdr[:1] + ["primary L_scan: 23701 alone for every module (the 23744 merge is in " + merged_out.name + ")"]
    else:
        hdr.append("23744: not merged (SS/out/lms_23744.npz absent); command: snake_scan.py lms --reduce scan "
                   "--npz23744 SS/out/lms_23744.npz")
    write_tsv(out, names, Lp, errp, np_, srcline, hdr, srcp, "source")
    return names, Lp, errp, np_, srcp, stats


def _reduce_with(f44, f01, names):
    """23744 without enough alpha events of its own: its reference normalisation uses the 23701 alpha medians."""
    f44 = dict(f44)
    f44["alpha_ref"] = np.concatenate([f["alpha_ref"] for f in f01])
    return reduce_group([f44])


def run_type_files():
    """{run type: [evio paths on disk]} from RUN_TYPE_OF (parsed from the benchmark source, not imported)."""
    src = BENCH.read_text()
    m = re.search(r"^RUN_TYPE_OF\s*=\s*(\{.*?\n\})", src, re.S | re.M)
    d = ast.literal_eval(m.group(1))
    out = {}
    for run, rt in sorted(d.items()):
        out.setdefault(rt, []).extend(sorted(glob.glob(str(EVIO_DIR / f"prad_{run:06d}.evio.*"))))
    return out


def reduce_prod(rts=None, out_dir=None):
    res = {}
    for rt, paths in run_type_files().items():
        if rts and rt not in rts:
            continue
        have = [(p, cache_path(p)) for p in paths if cache_path(p).exists()]
        if not have:
            continue
        files = [load(c) for _, c in have]
        names, L, err, n, info = reduce_group(files)
        runs = sorted({_run_split(p)[0] for p, _ in have})
        hdr = [f"run type {rt}: runs {runs}, {info['n_files']} of {len(paths)} local splits, "
               f"{info['n_lms']} LMS events, {info['n_alpha']} alpha events; "
               f"{info['n_files_group_alpha']} splits used the run type's alpha median (< {MIN_ALPHA} alpha events)"]
        write_tsv(Path(out_dir or OUT) / f"lms_prod_{rt}.tsv", names, L, err, n,
                  " ".join(os.path.basename(p) for p, _ in have) + "; map hycal_map.json (database)", hdr)
        res[rt] = (names, L, err, n, info, len(paths))
    return res


# ---------------------------------------------------------------------------------------------------------------
# dispatcher hooks
# ---------------------------------------------------------------------------------------------------------------

def add_parser(sub):
    p = sub.add_parser("lms", help="LMS pass: extraction per file, --reduce segments|scan|prod (D4, D8)")
    p.add_argument("evio", nargs="*")
    p.add_argument("--map", default=None, help="HyCal map; default: derived map for runs <= 23849, else database")
    p.add_argument("--cache-dir", default=str(NPZ_DIR))
    p.add_argument("--force", action="store_true")
    p.add_argument("--reduce", choices=("segments", "scan", "prod"))
    p.add_argument("--npz23744", default=str(OUT / "lms_23744.npz"))
    p.add_argument("--rt", nargs="*")
    return p


def run(a):
    if a.reduce == "segments":
        steps, rows, _ = reduce_segments()
        print(f"segments: {len(rows)} files, steps: {steps or 'none'}")
    elif a.reduce == "scan":
        _, L, err, n, src, st = reduce_scan(a.npz23744)
        print(f"scan: {int(np.isfinite(L).sum())} modules, median err {np.nanmedian(err):.4g}, "
              f"sources {dict(zip(*np.unique(src, return_counts=True)))}")
    elif a.reduce == "prod":
        for rt, (names, L, err, n, info, nl) in reduce_prod(a.rt).items():
            print(f"prod {rt}: {int(np.isfinite(L).sum())} modules, median err {np.nanmedian(err):.4g}, "
                  f"{info['n_lms']} LMS events in {info['n_files']}/{nl} splits")
    else:
        for e in a.evio:
            lms_pass(e, a.map if a.map is not None else default_map(e), cache_path(e, a.cache_dir), a.force)
    return 0
