#!/usr/bin/env python3
"""One-pulse-per-module reference templates, Wassim's algorithm (3.49 GeV H2 test).

Algorithm = Waveform_Wassim/best_ref_shape_fit.C (HyCalBestRef), via the existing
port in spline_template_extract.py (`select --wassim-ref`, `build --single-pulse`):
  candidate (macro evaluateCandidate, l.372-489): max sample > 2; peak time 150-170 ns;
    first sample > 2 before the peak with 20 samples before it (tier 2 if max|bg| <= 2,
    rms <= 1, roughness <= 1, else tier 1); >= 4 ascending samples above 2 up to the peak;
    >= 10 strictly falling samples after the peak; no local max > 2 and > 0.15*peak
    outside +-6 samples of the peak.
  ranking (betterCandidate, l.491-516): tier 2 over tier 1, then highest amplitude, then
    smaller roughness, then smaller bg rms; strict '<'/'>' so the first one wins ties.
  ONE pulse per module over all files; template = (W - bg_mean)/max, cubic interpolating
    spline through the 100 samples (TSpline3 analogue, l.872-943).  No pooling.
Literal readings / deviations (also in README.md):
  * template: unit peak = spline maximum (not max sample), onset shifted to 0, exp-tail
    pseudo-points to 416 ns (bank format needs T_END = 400 ns); macro: max sample = 1;
  * PbWO4 channels with 100 samples only (macro: every module_id);
  * pedestal = WaveAnalyzer pedestal (macro: ped_mean branch); samples in ADC (adc_to_mV=1);
  * extra non-macro veto: any sample >= 4095 (ADC clip) is skipped;
  * trigger gate physics (trigger_bits == 0x100); macro scans all events (max_events=-1 in its example);
  * modules without any candidate get the run-type average (bank needs every module);
  * the PDF omits the "h50 two-tau" overlay (requested); layout otherwise as write_pdf.
Usage: source setup_env.sh; python3 wassim_one_pulse_template.py [--workers 6]
"""
import argparse
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF = HERE / "output/templates_spline_best50_p1r3_20261104"
RUN_TYPE = "3p49GeV_H2"
OUT = HERE / "output/templates_wassim1_20261104" / RUN_TYPE


def evio_files():
    """The 3.49 GeV H2 split files of arm p1shape_r3 (from its select commands)."""
    pat = re.compile(r"select (\S+\.evio\.\d+) .*/p1shape_r3/%s/" % RUN_TYPE)
    return [m.group(1) for m in map(pat.search, (REF / "select_cmds.txt").read_text().splitlines()) if m]


def run_select(evio):
    name = re.sub(r".*prad_0*(\d+)\.evio\.(\d+)", r"\1_\2", evio)
    cmd = ["nice", "-n", "10", sys.executable, str(HERE / "spline_template_extract.py"), "select",
           evio, "--trigger", "physics", "--wassim-ref", "--out", str(OUT / "pulses" / f"{name}.npz")]
    with open(OUT / "logs" / f"select_{name}.log", "w") as log:
        return subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode


def build():
    """build --single-pulse in-process, with the two-tau overlay removed from the PDF pages."""
    import spline_template_extract as S
    orig = S._Pages.module

    def module_no_twotau(self, r, ctx):
        fig = orig(self, r, ctx)
        for ax in (self.a_main, self.a_tail):
            for ln in list(ax.lines):
                if ln.get_linestyle() == "--" and ln.get_color() == "k":   # the h50 two-tau curve
                    ln.remove()
        self.a_main.legend(fontsize=7, loc="upper right")
        return fig

    S._Pages.module = module_no_twotau
    return S.main(["build", "--single-pulse", "--pulses", str(OUT / "pulses" / "*.npz"),
                   "--run-type", RUN_TYPE, "--out", str(OUT / f"templates_{RUN_TYPE}")])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    t0 = time.monotonic()
    for d in ("pulses", "logs"):
        (OUT / d).mkdir(parents=True, exist_ok=True)
    files = evio_files()
    print(f"[select] {len(files)} files, {a.workers} workers", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        rc = list(ex.map(run_select, files))
    if any(rc):
        sys.exit(f"[ERROR] select failed for {sum(map(bool, rc))} file(s), see {OUT / 'logs'}")
    print(f"[select] done {time.monotonic() - t0:.0f} s", flush=True)
    sys.exit(build())
