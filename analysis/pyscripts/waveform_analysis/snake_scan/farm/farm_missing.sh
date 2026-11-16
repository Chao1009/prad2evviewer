#!/bin/bash
# farm_missing.sh CACHE_LIST OUTDIR   -- lines of CACHE_LIST (one /cache evio path per line) whose product OUTDIR/<run>/<file>.npz
# is absent: stdout is the list to resubmit (README step 4). Also removes stale partial files of dead jobs.
LIST="$1"; OUT="$2"
[ -s "$LIST" ] && [ -d "$OUT" ] || { echo "usage: $0 CACHE_LIST OUTDIR" >&2; exit 2; }
find "$OUT" -name '.*.partial.npz' -mmin +30 -delete 2>/dev/null
while IFS= read -r f; do
    b=$(basename "$f"); r=$((10#$(echo "$b" | sed -E 's/^prad_0*([0-9]+)\..*/\1/')))
    [ -s "$OUT/$r/$b.npz" ] || echo "$f"
done < "$LIST"
