#!/bin/bash
# jcache.sh [LIST] [--check|--dry-run]   -- README step 2b (run on ifarm, after make_filelist.sh).
#   jcache.sh jcache_list.txt            request tape -> /cache for every /mss path of the list (jcache get, batches of $BATCH)
#   jcache.sh jcache_list.txt --check    count how many files are already on /cache (exit 1 if any is missing); jcache get is
#                                        asynchronous: re-run --check until it reports 0 missing before submitting
#   jcache.sh jcache_list.txt --dry-run  print the commands only
# Cache location: /mss/clas12/rg-o/data/X  ->  /cache/clas12/rg-o/data/X  (the usual JLab mapping; override with CACHE_ROOT/MSS_ROOT).
# The cache is a scratch area with an automatic cleanup: submit soon after the files arrive (jcache pin is optional, see
# `jcache help` on ifarm for its option syntax; not used here).
set -u
LIST="${1:-jcache_list.txt}"; MODE="${2:-get}"
MSS_ROOT="${MSS_ROOT:-/mss/clas12/rg-o/data}"; CACHE_ROOT="${CACHE_ROOT:-/cache/clas12/rg-o/data}"
BATCH="${BATCH:-100}"
[ -s "$LIST" ] || { echo "list $LIST missing or empty (run make_filelist.sh)" >&2; exit 2; }
case "$MODE" in
  get|--get)
    echo "requesting $(wc -l < "$LIST") files from tape"
    split -l "$BATCH" "$LIST" "${TMPDIR:-/tmp}/jcache_batch_$$_"
    for b in "${TMPDIR:-/tmp}/jcache_batch_$$_"*; do jcache get $(cat "$b"); rm -f "$b"; done
    echo "done; follow with: $0 $LIST --check";;
  --dry-run)
    split -l "$BATCH" "$LIST" "${TMPDIR:-/tmp}/jcache_batch_$$_"
    for b in "${TMPDIR:-/tmp}/jcache_batch_$$_"*; do echo "jcache get $(tr '\n' ' ' < "$b" | cut -c1-200)..."; rm -f "$b"; done;;
  --check|check)
    miss=0; have=0
    while IFS= read -r f; do
        c="${f/#$MSS_ROOT/$CACHE_ROOT}"
        if [ -e "$c" ]; then have=$((have+1)); else miss=$((miss+1)); [ "$miss" -le 10 ] && echo "missing: $c"; fi
    done < "$LIST"
    echo "on /cache: $have, missing: $miss"
    [ "$miss" -eq 0 ];;
  *) echo "usage: $0 LIST [--check|--dry-run]" >&2; exit 2;;
esac
