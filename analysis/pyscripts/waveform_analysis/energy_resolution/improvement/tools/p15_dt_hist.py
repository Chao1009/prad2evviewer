"""P15 figure: seed-minus-member time dt = t_seed - t_member for non-seed members of legacy (ungated) clusters, raw vs offset-corrected.
Source: crosstalk_study/cache/<run type>/*_legacy.npz (first NSPLIT splits per cache dir; Ta pools its three run dirs). No e-p cut.
calibrated = raw - offset_ns (time_calib.json, resolved rows; unresolved and lead glass 0). Usage: python p15_dt_hist.py"""
import glob, json, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PY = HERE.parents[3]                      # analysis/pyscripts
CACHE = PY / "crosstalk_study" / "cache"
DB = Path(os.environ.get("PRAD2_DATABASE_DIR", PY.parents[1] / "database"))
OUT = HERE.parent / "figures" / "P15_seed_member_dt"
OUT.parent.mkdir(exist_ok=True)
NSPLIT = 3
PANELS = [("0p73GeV_H2", ["0p73GeV_H2"], 12, "0.73 GeV H2"), ("2p24GeV_H2", ["2p24GeV_H2"], 16, "2.24 GeV H2"),
          ("2p24GeV_Ta_X17", ["2p24GeV_Ta_X17_25800", "2p24GeV_Ta_X17_26138", "2p24GeV_Ta_X17_26230"], 4, "2.24 GeV Ta X17"),
          ("3p49GeV_H2", ["3p49GeV_H2"], 12, "3.49 GeV H2")]
BINS = np.arange(-20, 30.01, 0.5)
C_RAW, C_CAL, C_GATE = "#7a7a7a", "#1f6fb4", "#c2410c"


def offsets():
    j = json.load(open(DB / "hycal_time_offsets" / "time_calib.json"))["modules"]
    return {e["name"]: float(e["offset_ns"]) for e in j if e.get("resolved")}


def rsig(x):
    lo, hi = np.percentile(x, [15.865, 84.135])
    return 0.5 * (hi - lo)


def load(dirs, off):
    raw, cal = [], []
    for d in dirs:
        for f in sorted(glob.glob(str(CACHE / d / "*_legacy.npz")))[:NSPLIT]:
            z = np.load(f)
            names = [str(n) for n in z["module_names"]]
            o = np.array([off.get(n, 0.0) for n in names])
            mi, si = z["mod_i"].astype(np.int64), z["seed_mod_i"].astype(np.int64)
            ts, t = z["seed_t"].astype(np.float64), z["t"].astype(np.float64)
            raw.append(ts - t)
            cal.append((ts - o[si]) - (t - o[mi]))
    return np.concatenate(raw), np.concatenate(cal)


def main():
    off = offsets()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    for ax, (key, dirs, gate, label) in zip(axs.flat, PANELS):
        raw, cal = load(dirs, off)
        sr, sc = rsig(raw), rsig(cal)
        print(f"{key}: N={raw.size} robust sigma raw {sr:.2f} cal {sc:.2f} ns; frac dt>8: raw {np.mean(raw>8):.4f} cal {np.mean(cal>8):.4f}; "
              f"dt<-8: raw {np.mean(raw<-8):.4f} cal {np.mean(cal<-8):.4f}; median raw {np.median(raw):.2f} cal {np.median(cal):.2f}")
        ax.hist(raw, BINS, histtype="step", color=C_RAW, lw=1.4, label=f"raw, robust $\\sigma$ = {sr:.2f} ns")
        ax.hist(cal, BINS, histtype="step", color=C_CAL, lw=1.6, label=f"offset-corrected, robust $\\sigma$ = {sc:.2f} ns")
        ax.axvline(gate, color=C_GATE, ls="--", lw=1.2, label=f"P15 best gate {gate} ns")
        ax.set_yscale("log")
        ax.set_xlim(-20, 30)
        ax.set_title(f"{label}  ({raw.size/1e6:.2f} M members)", fontsize=11)
        ax.grid(alpha=0.2, lw=0.5)
        ax.legend(frameon=False, fontsize=9, loc="upper left")
        top = max(np.histogram(raw, BINS)[0].max(), np.histogram(cal, BINS)[0].max())
        ax.set_ylim(0.5, top * 60)
    for ax in axs[1]:
        ax.set_xlabel(r"$\Delta t = t_{seed} - t_{member}$ (ns)")
    for ax in axs[:, 0]:
        ax.set_ylabel("members / 0.5 ns")
    fig.suptitle(r"P15: seed-to-member time difference, legacy (ungated) clusters, non-seed members", fontsize=13)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}.{ext}", dpi=150)


if __name__ == "__main__":
    main()
