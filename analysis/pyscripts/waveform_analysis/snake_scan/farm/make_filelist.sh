#!/bin/bash
# make_filelist.sh [OUT]   -- run on ifarm (needs the /mss tape-stub listing). README step 2a.
# Builds the EXACT jcache list from what is on tape, not from a guessed file count:
#   scan runs 23699-23723: every STEP-th split, i.e. splits whose 5-digit number NNNNN satisfies NNNNN % STEP == OFFSET
#                          (STEP=3, OFFSET=0: 00000, 00003, 00006, ...); a split missing from tape is simply absent
#   run 23744 (cosmic)   : all splits
# Output: OUT (default jcache_list.txt, one /mss path per line) and OUT with .txt replaced by _cache.txt (same files under /cache).
# Settings (env): MSS_ROOT (default /mss/clas12/rg-o/data), SCAN_RUNS, COSMIC_RUN, STEP, OFFSET.
set -u
OUT="${1:-jcache_list.txt}"
MSS_ROOT="${MSS_ROOT:-/mss/clas12/rg-o/data}"
SCAN_RUNS="${SCAN_RUNS:-$(seq 23699 23723 | tr '\n' ' ')}"
COSMIC_RUN="${COSMIC_RUN:-23744}"
STEP="${STEP:-3}"; OFFSET="${OFFSET:-0}"
CACHE_ROOT="${CACHE_ROOT:-/cache/clas12/rg-o/data}"
OUTC="${OUT%.txt}_cache.txt"
: > "$OUT"; : > "$OUTC"
nbad=0
printf "%-8s %8s %9s\n" run on_tape selected
list_run() {   # run  mode(step|all)
    local run="$1" mode="$2" dir="$MSS_ROOT/prad_0$1" n=0 sel=0 f nn
    if [ ! -d "$dir" ]; then echo "$run: directory $dir not found" >&2; nbad=$((nbad+1)); return; fi
    while IFS= read -r f; do
        n=$((n+1))
        nn=$((10#${f##*.}))
        if [ "$mode" = all ] || [ $((nn % STEP)) -eq "$OFFSET" ]; then
            echo "$dir/$f" >> "$OUT"; echo "$CACHE_ROOT/prad_0$run/$f" >> "$OUTC"; sel=$((sel+1))
        fi
    done < <(ls "$dir" | grep -E "^prad_0${run}\.evio\.[0-9]{5}$" | sort)
    printf "%-8s %8d %9d\n" "$run" "$n" "$sel"
}
for r in $SCAN_RUNS; do list_run "$r" step; done
list_run "$COSMIC_RUN" all
echo "wrote $(wc -l < "$OUT" | tr -d " ") paths to $OUT and $OUTC ($nbad run directories missing)"
[ "$nbad" -eq 0 ]
