"""
EVIO file helpers shared by the GUIs (Qt-free).

* Fetching raw files from the DAQ host: default locations, local-file
  listing, disk-space estimate and the scp copy script.
* Record access through prad2py: open a file, walk its physics records,
  and seek to a record by index.  prad2py is imported on first use, so a
  GUI can start and report a missing prad2py itself.
"""

from __future__ import annotations

import glob
import itertools
import os
import re
import shutil
import subprocess
from typing import Callable, Iterator, List, Optional, Tuple

REMOTE_HOST = "clondaq2"
REMOTE_DATA_BASE = "/data/stage2"
LOCAL_DATA_BASE = "/data/evio/data"
# Conservative per-file size for disk-space estimates when an actual remote
# `ls -l` listing isn't available.  Each evio file is ~2 GB; bump slightly
# for safety so we don't run out mid-copy.
EVIO_BYTES_PER_FILE_EST = int(2.1 * 1024 ** 3)

_EVIO_INDEX_RE = re.compile(r'\.evio\.(\d+)$')


# ---- Remote fetch ----------------------------------------------------------

def evio_index(name: str) -> Optional[int]:
    """File number N of a ``*.evio.N`` name, or None."""
    m = _EVIO_INDEX_RE.search(name)
    return int(m.group(1)) if m else None


def local_evio_in_range(run_dir: str, f_start: int, f_end: int,
                        pattern: str = "*.evio.*") -> List[str]:
    """Sorted basenames of the files in ``run_dir`` matching ``pattern``
    whose file number lies in [f_start, f_end] (empty if the directory does
    not exist)."""
    out: List[str] = []
    for path in sorted(glob.glob(os.path.join(run_dir, pattern))):
        name = os.path.basename(path)
        n = evio_index(name)
        if n is not None and f_start <= n <= f_end:
            out.append(name)
    return out


def free_bytes(path: str) -> int:
    """Free space on the filesystem holding ``path``, or its nearest
    existing ancestor when ``path`` does not exist yet."""
    check_path = path
    while check_path and not os.path.exists(check_path):
        check_path = os.path.dirname(check_path)
    return shutil.disk_usage(check_path or "/").free


def fmt_bytes(b: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if b < 1024:
            return f"{b:.1f} {unit}"
        b /= 1024
    return f"{b:.1f} PB"


def check_disk_space(remote_host: str, remote_run_dir: str,
                     local_base: str, f_start: int, f_end: int,
                     local_run_dir: Optional[str] = None) -> Tuple[int, int]:
    """Return (needed_bytes, free_bytes) for evio files [f_start, f_end].

    SSHes to remote_host and sums the sizes of files whose .evio.NNN suffix
    falls within [f_start, f_end].  Files that already exist in local_run_dir
    are excluded from the calculation.  If the remote listing has no file in
    the range (remote directory missing or still empty), falls back to a
    conservative ~2 GB-per-file estimate for the file numbers not present
    locally.  Free space is measured on the filesystem that contains
    local_base (or its nearest existing ancestor).
    Raises RuntimeError if the SSH call itself fails (exit 255).
    """
    result = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=10",
         remote_host, f"ls -l {remote_run_dir}/ 2>/dev/null"],
        capture_output=True, text=True, timeout=30,
    )
    if result.returncode == 255:
        raise RuntimeError(result.stderr.strip() or "SSH connection failed")

    needed = 0
    listed = 0
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) < 9:
            continue
        fname = parts[-1]
        n = evio_index(fname)
        if n is None or not f_start <= n <= f_end:
            continue
        listed += 1
        if local_run_dir and os.path.isfile(os.path.join(local_run_dir, fname)):
            continue  # already downloaded
        try:
            needed += int(parts[4])
        except ValueError:
            needed += EVIO_BYTES_PER_FILE_EST

    if listed == 0:
        present = set()
        if local_run_dir:
            present = {evio_index(f) for f in
                       local_evio_in_range(local_run_dir, f_start, f_end)}
        missing = sum(1 for n in range(f_start, f_end + 1)
                      if n not in present)
        needed = missing * EVIO_BYTES_PER_FILE_EST

    return needed, free_bytes(local_base)


def scp_bash(remote_host: str, remote_run_dir: str, local_run_dir: str,
             f_start: int, f_end: int) -> str:
    """Bash script that lists ``remote_run_dir`` on ``remote_host`` and
    scps the evio files numbered [f_start, f_end] that are not yet in
    ``local_run_dir``."""
    return (
        f"mkdir -p {local_run_dir}\n"
        f"echo 'Local directory: {local_run_dir}'\n"
        f"echo 'Listing remote files...'\n"
        f"ALL_FILES=$(ssh {remote_host} 'ls {remote_run_dir}/' 2>/dev/null | sort)\n"
        f"COPIED=0\n"
        f"ALREADY=0\n"
        f"while IFS= read -r f; do\n"
        f"    NUM=$(echo \"$f\" | grep -oP '\\.evio\\.\\K[0-9]+')\n"
        f"    [ -z \"$NUM\" ] && continue\n"
        f"    N=$((10#$NUM))\n"
        f"    if [ \"$N\" -lt {f_start} ] || [ \"$N\" -gt {f_end} ]; then continue; fi\n"
        f"    if [ -f \"{local_run_dir}/$f\" ]; then\n"
        f"        echo \"  Already exists: $f (skipping)\"\n"
        f"        ALREADY=$((ALREADY+1))\n"
        f"    else\n"
        f"        echo \"  Copying $f\"\n"
        f"        scp {remote_host}:{remote_run_dir}/$f {local_run_dir}/\n"
        f"        COPIED=$((COPIED+1))\n"
        f"    fi\n"
        f"done <<< \"$ALL_FILES\"\n"
        f"echo \"Done. Copied $COPIED file(s), $ALREADY already present.\"\n"
    )


# ---- Record access (prad2py) -----------------------------------------------

def _dec():
    from prad2py import dec
    return dec


def _daq_config(daq_cfg):
    """DaqConfig from a DaqConfig, a daq_config.json path, or None/"" for
    the installed default."""
    if daq_cfg is None or isinstance(daq_cfg, (str, os.PathLike)):
        return _dec().load_daq_config(os.fspath(daq_cfg) if daq_cfg else "")
    return daq_cfg


def open_evio(path, daq_cfg=None):
    """Open ``path`` with EvChannel::OpenAuto (random access when the file
    supports it, else sequential).

    ``daq_cfg`` is a DaqConfig, a daq_config.json path, or None for the
    installed default.  Returns ``(channel, is_random_access)``; raises
    RuntimeError if the file cannot be opened.
    """
    dec = _dec()
    ch = dec.EvChannel()
    ch.set_config(_daq_config(daq_cfg))
    st = ch.open_auto(str(path))
    if st != dec.Status.success:
        raise RuntimeError(f"cannot open {path}: {st}")
    return ch, bool(ch.is_random_access())


def iter_physics_records(ch, is_ra: bool,
                         cancel: Optional[Callable[[], bool]] = None,
                         on_record: Optional[Callable[[int], None]] = None,
                         ) -> Iterator[int]:
    """Walk every record of an open channel and yield the 0-based index of
    each Physics record, with ``ch`` holding it already scanned.

    Random-access mode visits every index and skips records that fail to
    read; sequential mode stops at the first failed read (EOF).  ``cancel``
    is polled before each record; ``on_record(idx)`` runs after each record,
    physics or not, once the consumer is done with it.
    """
    dec = _dec()
    indices = (range(ch.get_random_access_event_count()) if is_ra
               else itertools.count())
    for idx in indices:
        if cancel is not None and cancel():
            return
        if is_ra:
            ok = ch.read_event_by_index(idx) == dec.Status.success
        elif ch.read() == dec.Status.success:
            ok = True
        else:
            return
        if ok and ch.scan() and ch.get_event_type() == dec.EventType.Physics:
            yield idx
        if on_record is not None:
            on_record(idx)


class EvioCursor:
    """Positioned reader over one file: ``seek(rec)`` loads record ``rec``
    (0-based, as yielded by iter_physics_records) into ``ch``, unscanned.

    Random access jumps directly.  Sequential mode walks forward, reopens
    the file for a backward seek, and leaves the buffer alone when ``rec``
    is the record already loaded, so sub-events of one record can be
    visited back to back.  Only use ``ch`` on the loaded record (scan,
    select_event, decode...): reading it directly desynchronises a
    sequential cursor.
    """

    def __init__(self, path, daq_cfg=None):
        self.path = str(path)
        self._cfg = _daq_config(daq_cfg)
        self.ch, self.is_ra = open_evio(self.path, self._cfg)
        self._pos = -1

    def seek(self, rec: int) -> None:
        """Raises RuntimeError if the record cannot be read."""
        dec = _dec()
        if self.is_ra:
            st = self.ch.read_event_by_index(rec)
            if st != dec.Status.success:
                raise RuntimeError(f"read_event_by_index({rec}) failed: {st}")
            return
        if self._pos > rec:
            ch, _ = open_evio(self.path, self._cfg)
            self.ch.close()
            self.ch, self._pos = ch, -1
        while self._pos < rec:
            if self.ch.read() != dec.Status.success:
                raise RuntimeError(f"EOF before reaching record {rec}")
            self._pos += 1

    def close(self) -> None:
        self.ch.close()
