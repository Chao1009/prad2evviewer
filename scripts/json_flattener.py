#!/usr/bin/env python3
"""
JSON Flattener
==============
Flatten a JSON document into a flat (key, value) table.

Nested objects become dotted paths; arrays use numeric indices:

    {"a": {"b": 1}, "c": [10, 20]}
        →
        a.b   1
        c.0   10
        c.1   20

Depth limit
-----------
``--max-depth N`` caps recursion. Depth = number of key segments in
the emitted path:

    {"a": {"b": {"c": 1}}}     max-depth=1  →  a       {"b":{"c":1}}
                               max-depth=2  →  a.b     {"c":1}
                               max-depth=3  →  a.b.c   1
                               max-depth=-1 →  a.b.c   1     (default, unlimited)

When a container is reached at the cap it is kept as a compact JSON
string and one warning per such value is printed to stderr. Scalars
at or past the cap are emitted normally — only containers trigger
truncation.

Usage
-----
    python scripts/json_flattener.py input.json
    python scripts/json_flattener.py input.json --max-depth 2
    python scripts/json_flattener.py input.json --separator /
    python scripts/json_flattener.py input.json --csv > out.csv
    cat input.json | python scripts/json_flattener.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import Any


def flatten(obj: Any, max_depth: int = -1, separator: str = "."):
    """Return (rows, truncated_paths) for the flattened JSON."""
    rows: list[tuple[str, Any]] = []
    truncated: list[str] = []

    def walk(node: Any, path: str, depth: int) -> None:
        is_container = isinstance(node, (dict, list))
        if max_depth >= 0 and depth >= max_depth and is_container:
            truncated.append(path or "<root>")
            rows.append((path, json.dumps(node, separators=(",", ":"),
                                          ensure_ascii=False)))
            return

        if isinstance(node, dict):
            if not node:
                rows.append((path, "{}"))
                return
            for k, v in node.items():
                child = f"{path}{separator}{k}" if path else str(k)
                walk(v, child, depth + 1)
        elif isinstance(node, list):
            if not node:
                rows.append((path, "[]"))
                return
            for i, v in enumerate(node):
                child = f"{path}{separator}{i}" if path else str(i)
                walk(v, child, depth + 1)
        else:
            rows.append((path, node))

    walk(obj, "", 0)
    return rows, truncated


def value_to_str(v: Any) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def _tsv_escape(s: str) -> str:
    return s.replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Flatten a JSON document into a (key, value) table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  json_flattener.py database/config.json\n"
            "  json_flattener.py runinfo.json --max-depth 2 --csv > runs.csv\n"
            "  curl -s http://host/state.json | json_flattener.py -\n"
        ),
    )
    p.add_argument("input", nargs="?", default="-",
                   help="JSON file path, or '-' for stdin (default).")
    p.add_argument("--max-depth", type=int, default=-1,
                   help="Max nesting depth; deeper containers become "
                        "JSON-string leaves with a warning. -1 = unlimited "
                        "(default).")
    p.add_argument("--separator", default=".",
                   help="Path-component joiner (default '.').")
    p.add_argument("--csv", action="store_true",
                   help="Emit CSV instead of TSV.")
    p.add_argument("--no-header", action="store_true",
                   help="Skip the 'key/value' header row.")
    args = p.parse_args(argv)

    try:
        if args.input == "-":
            data = json.load(sys.stdin)
        else:
            with open(args.input, "r", encoding="utf-8") as f:
                data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        sys.stderr.write(f"[ERROR] cannot read JSON from "
                         f"{args.input!r}: {e}\n")
        return 1

    rows, truncated = flatten(data, max_depth=args.max_depth,
                              separator=args.separator)

    if args.csv:
        w = csv.writer(sys.stdout, lineterminator="\n")
        if not args.no_header:
            w.writerow(["key", "value"])
        for k, v in rows:
            w.writerow([k, value_to_str(v)])
    else:
        if not args.no_header:
            sys.stdout.write("key\tvalue\n")
        for k, v in rows:
            sys.stdout.write(f"{k}\t{_tsv_escape(value_to_str(v))}\n")

    if truncated:
        sys.stderr.write(
            f"[WARN] json_flattener: {len(truncated)} value(s) kept as "
            f"JSON string at max-depth={args.max_depth}:\n"
        )
        for path in truncated:
            sys.stderr.write(f"  - {path}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
