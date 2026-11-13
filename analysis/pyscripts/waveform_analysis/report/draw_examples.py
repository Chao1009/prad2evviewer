"""3 smooth + 3 bumpy example tails from the pooled scan reservoirs (smallest random keys)."""
import glob, numpy as np
import matplotlib.pyplot as plt
# plot style: large bold axis labels, large legends
plt.rcParams.update({
    "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
    "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
files = sorted(glob.glob("output/pulse_shape_study/scan/run_*.npz"))
zs = [np.load(f, allow_pickle=False) for f in files]
names = zs[0]["mod_name"]
keys = ["midx", "event_number", "amp", "peak", "cls", "dip", "bump_adc", "bump", "gap", "samples", "ped", "key"]
ex = {k: np.concatenate([z["ex_" + k] for z in zs]) for k in keys}
runs = np.concatenate([[int(z["run"])] * z["ex_key"].size for z in zs])
# modules.tsv top list
rows = [l.rstrip("\n").split("\t") for l in open("output/pulse_shape_study/modules.tsv")][1:]
top = sorted([(r[0], int(r[7]), int(r[8]), float(r[9])) for r in rows if int(r[7]) >= 100], key=lambda x: -x[3])[:6]
print("top n>=100:", top)
print("W554:", [r for r in rows if r[0] == "W554"])
fig, axs = plt.subplots(2, 3, figsize=(18, 10))
for row, (cls, label) in enumerate(((0, "smooth"), (2, "bumpy"))):
    ix = np.nonzero(ex["cls"] == cls)[0]
    ix = ix[np.argsort(ex["key"][ix])[:3]]           # uniform sample of the class population
    for ax, i in zip(axs[row], ix):
        w = ex["samples"][i].astype(float) - ex["ped"][i]
        pk, dip, gap = int(ex["peak"][i]), int(ex["dip"][i]), int(ex["gap"][i])
        lo, hi = max(pk - 4, 0), min(pk + 12, w.size - 1)
        t = np.arange(lo, hi + 1)
        ax.plot(t, w[lo:hi + 1], "-", color="#333", lw=1.2, marker=".", ms=4)
        ax.axvspan(pk + 0.5, pk + 8.5, color="#999", alpha=0.12, lw=0, label="tail window")
        ax.plot(pk, w[pk], "o", color="tab:blue", ms=7, label="peak")
        if dip >= 0:
            ax.plot(dip, w[dip], "v", color="tab:red", ms=8, label="dip")
            ax.plot(pk + gap, w[pk + gap], "^", color="tab:orange", ms=8, label="secondary max")
        ax.set_ylim(0, w[pk] * 1.12)
        ax.grid(True, color="#ddd", lw=0.5)
        ttl = f"{names[ex['midx'][i]]}  run {runs[i]}  ev {ex['event_number'][i]}\n{label}: A = {ex['amp'][i]:.0f}"
        if cls == 2:
            ttl += f", bump = {ex['bump_adc'][i]:.0f} ADC ({100 * ex['bump'][i]:.1f} % A), gap = {gap}"
        ax.set_title(ttl, fontsize=11)
        print(label, names[ex['midx'][i]], runs[i], ex['event_number'][i], f"A={ex['amp'][i]:.0f}", f"bump={ex['bump_adc'][i]:.1f}", f"{100*ex['bump'][i]:.1f}%", "dip-peak", dip - pk if dip >= 0 else None, "gap", gap)
for ax in axs[-1]:
    ax.set_xlabel("sample (4 ns)")
for ax in axs[:, 0]:
    ax.set_ylabel("ADC above pedestal")
h, l = axs[1, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Example tails, A >= 200 ADC: 3 smooth (top), 3 bumpy (bottom); shaded = 8-sample tail window", fontsize=14)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig("waveform_analysis/report/figures/examples.png", dpi=150)
