#!/usr/bin/env python3
"""Peaks-only raw ROOT writer: MF pulse energies in the WaveAnalyzer peak slots.

Reads a standard raw ROOT file that still holds waveforms (prad2ana_replay_rawdata,
default mode), runs the spline-template matched filter (mf_pulse.py, R5nr settings)
on its PbWO4 channels and writes a raw file with the branches that
prad2ana_physics_calib reads, where hycal.peak_integral is E_u = A * I_int * f(A) * kappa
(kappa = 1) in ADC*samples, hycal.peak_height = A and hycal.peak_time = t0 + t50 (raw ns).
The unchanged physics_calib then calibrates on MF energies (plans/plan_mf_stdcal_20261102.md).

Fallbacks follow B (benchmark_hycal_timing_test_20261021.py): non-PbWO4 channels, channels
with an unexpected sample count and clipped PbWO4 channels flagged wa_fb take the
WaveAnalyzer peaks (integral as is, no f(A)).  peak_quality bit 7 marks an MF peak.

Only sum-trigger events with hycal.nch <= 100 are written (physics_calib skips the rest).
The output keeps the input basename.  uproot cannot write the [nch][8] leaves in a form
ROOT reads (it stores them as jagged), so the peaks go through a flat uproot tree and a
small ACLiC macro (root must be on PATH) that fills the real branches.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

import _common as C   # noqa: F401  (sets up prad2py import message)
import mf_pulse
from prad2py import det

TBIT_SUM = 1 << 8
NCH_MAX = 100
MAX_PEAKS = 8
MF_BIT = 0x80

MACRO = r'''
#include "TFile.h"
#include "TTree.h"
void fill_peaks(const char *fin, const char *fout)
{
    const int MC = 1024;
    TFile fi(fin);
    TTree *ti = (TTree *)fi.Get("flat");
    Int_t event_num, nch, n8;
    UInt_t trigger_bits;
    UChar_t trigger_type;
    static UShort_t mid[MC];
    static UChar_t mtype[MC], npk[MC], pq[MC * 8];
    static Float_t gf[MC], ph[MC * 8], pt[MC * 8], pi[MC * 8];
    ti->SetBranchAddress("event_num", &event_num);
    ti->SetBranchAddress("trigger_type", &trigger_type);
    ti->SetBranchAddress("trigger_bits", &trigger_bits);
    ti->SetBranchAddress("nch", &nch);
    ti->SetBranchAddress("n8", &n8);
    ti->SetBranchAddress("module_id", mid);
    ti->SetBranchAddress("module_type", mtype);
    ti->SetBranchAddress("gain_factor", gf);
    ti->SetBranchAddress("npeaks", npk);
    ti->SetBranchAddress("pk_h", ph);
    ti->SetBranchAddress("pk_t", pt);
    ti->SetBranchAddress("pk_i", pi);
    ti->SetBranchAddress("pk_q", pq);

    TFile fo(fout, "RECREATE");
    TTree *to = new TTree("events", "events");
    Int_t o_event_num, o_nch;
    UInt_t o_tbits;
    UChar_t o_ttype;
    static UShort_t o_mid[MC];
    static UChar_t o_mtype[MC], o_npk[MC], o_pq[MC][8];
    static Float_t o_gf[MC], o_ph[MC][8], o_pt[MC][8], o_pi[MC][8];
    to->Branch("event_num", &o_event_num, "event_num/I");
    to->Branch("trigger_type", &o_ttype, "trigger_type/b");
    to->Branch("trigger_bits", &o_tbits, "trigger_bits/i");
    to->Branch("hycal.nch", &o_nch, "hycal.nch/I");
    to->Branch("hycal.module_id", o_mid, "hycal.module_id[hycal.nch]/s");
    to->Branch("hycal.module_type", o_mtype, "hycal.module_type[hycal.nch]/b");
    to->Branch("hycal.gain_factor", o_gf, "hycal.gain_factor[hycal.nch]/F");
    to->Branch("hycal.npeaks", o_npk, "hycal.npeaks[hycal.nch]/b");
    to->Branch("hycal.peak_height", o_ph, "hycal.peak_height[hycal.nch][8]/F");
    to->Branch("hycal.peak_time", o_pt, "hycal.peak_time[hycal.nch][8]/F");
    to->Branch("hycal.peak_integral", o_pi, "hycal.peak_integral[hycal.nch][8]/F");
    to->Branch("hycal.peak_quality", o_pq, "hycal.peak_quality[hycal.nch][8]/b");
    Long64_t n = ti->GetEntries();
    for (Long64_t e = 0; e < n; ++e) {
        ti->GetEntry(e);
        if (nch > MC) continue;
        o_event_num = event_num; o_ttype = trigger_type; o_tbits = trigger_bits;
        o_nch = nch;
        for (int j = 0; j < nch; ++j) {
            o_mid[j] = mid[j]; o_mtype[j] = mtype[j]; o_gf[j] = gf[j]; o_npk[j] = npk[j];
            for (int k = 0; k < 8; ++k) {
                o_ph[j][k] = ph[j * 8 + k]; o_pt[j][k] = pt[j * 8 + k];
                o_pi[j][k] = pi[j * 8 + k]; o_pq[j][k] = pq[j * 8 + k];
            }
        }
        to->Fill();
    }
    to->Write();
    fo.Close();
}
'''


def sha1_file(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def wa_peaks(wave_ana, samples, ns):
    s = np.ascontiguousarray(samples[:ns], dtype=np.uint16)
    _, _, pks = wave_ana.analyze(s)
    return [(pk.time, pk.height, pk.integral) for pk in pks]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input", help="raw ROOT file with hycal.samples")
    ap.add_argument("--mf-templates", required=True, help="spline_templates_v1 .npz (R5nr bank)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config JSON (default: PipelineBuilder default)")
    ap.add_argument("--mf-r-model", type=float, default=0.3)
    ap.add_argument("--mf-min-amp", type=float, default=5.0)
    ap.add_argument("--mf-kmax", type=int, default=4)
    ap.add_argument("--mf-nsigma", type=float, default=5.0)
    ap.add_argument("--max-events", type=int, default=0, help="stop after N input entries")
    ap.add_argument("--chunk", type=int, default=4000)
    ap.add_argument("--diag", action="store_true",
                    help="also write <stem>_ratio.tsv (per-module E_u / WA integral)")
    ap.add_argument("--diag-events", type=int, default=30000,
                    help="--diag: use the first N input entries only")
    args = ap.parse_args()

    import uproot
    import awkward as ak

    t_start = time.time()
    in_path = os.path.abspath(args.input)
    base = os.path.basename(in_path)
    stem = base[:-5] if base.endswith(".root") else base
    os.makedirs(args.out, exist_ok=True)
    out_path = os.path.join(os.path.abspath(args.out), base)

    m = re.search(r"prad_0*(\d+)", base)
    run = int(m.group(1)) if m else -1
    b = det.PipelineBuilder()
    if args.daq_config:
        b.set_daq_config(args.daq_config)
    if run > 0:
        b.set_run_number(run)
    core = b.build()
    hycal = core.hycal
    from prad2py import dec
    wave_ana = dec.WaveAnalyzer(core.daq_cfg.wave_cfg)
    wcfg = wave_ana.cfg
    mods = [hycal.module(i) for i in range(hycal.module_count())]
    names = [mm.name for mm in mods]
    pwo4 = [bool(mm.is_pwo4()) for mm in mods]
    idx_of_id = np.full(65536, -1, dtype=np.int64)
    for mm in mods:
        idx_of_id[mm.id] = mm.index
    pw_arr = np.array(pwo4, dtype=bool)

    bank, info = mf_pulse.build_bank(
        names, pwo4, args.mf_templates,
        dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)))
    rec = mf_pulse.MFReconstructor(
        bank, nsigma=args.mf_nsigma, min_amp=args.mf_min_amp, kmax=args.mf_kmax,
        r_model=args.mf_r_model, ped_nsamples=int(wcfg.ped_nsamples))
    overflow = int(wcfg.overflow)

    st = dict(n_events=0, n_in=0, n_skip_trig=0, n_skip_nch=0, n_chan=0, n_mf=0,
              n_wa_fb=0, n_wa_fb_peaks=0, n_odd=0, n_odd_peaks=0, n_lg=0, n_lg_peaks=0,
              n_peak_chan=0, n_dropped=0, n_mf_pulses=0, n_sat_fallback=0)
    ratio = {}   # module name -> (list of ratio, list of dt)

    tmpd = tempfile.mkdtemp(prefix="sup-mfraw-", dir=os.path.abspath(args.out))
    flat_path = os.path.join(tmpd, "flat.root")
    flat = uproot.recreate(flat_path)
    types = dict(event_num=np.int32, trigger_type=np.uint8, trigger_bits=np.uint32,
                 nch=np.int32, n8=np.int32, module_id="var * uint16",
                 module_type="var * uint8", gain_factor="var * float32",
                 npeaks="var * uint8", pk_h="var * float32", pk_t="var * float32",
                 pk_i="var * float32", pk_q="var * uint8")
    flat.mktree("flat", types,
                counter_name=lambda c: "n8" if c.startswith("pk_") else "nch")

    tin = uproot.open(in_path + ":events")
    n_entries = tin.num_entries
    if args.max_events:
        n_entries = min(n_entries, args.max_events)
    keys = ["event_num", "trigger_type", "trigger_bits", "hycal.nch", "hycal.module_id",
            "hycal.module_type", "hycal.nsamples", "hycal.samples", "hycal.gain_factor"]
    NS = None
    for start in range(0, n_entries, args.chunk):
        stop = min(start + args.chunk, n_entries)
        a = tin.arrays(keys, entry_start=start, entry_stop=stop, library="ak")
        st["n_in"] += stop - start
        tb = ak.to_numpy(a["trigger_bits"])
        nchv = ak.to_numpy(a["hycal.nch"])
        sum_ok = (tb & TBIT_SUM) > 0
        sel = sum_ok & (nchv <= NCH_MAX)
        st["n_skip_trig"] += int((~sum_ok).sum())
        st["n_skip_nch"] += int((sum_ok & ~sel).sum())
        if not sel.any():
            continue
        a = a[sel]
        nch = ak.to_numpy(a["hycal.nch"]).astype(np.int32)
        ne = nch.size
        evn = ak.to_numpy(a["event_num"]).astype(np.int32)
        ttype = ak.to_numpy(a["trigger_type"]).astype(np.uint8)
        tbits = ak.to_numpy(a["trigger_bits"]).astype(np.uint32)
        mid_id = ak.to_numpy(ak.flatten(a["hycal.module_id"])).astype(np.int64)
        mtype = ak.to_numpy(ak.flatten(a["hycal.module_type"])).astype(np.uint8)
        nsm = ak.to_numpy(ak.flatten(a["hycal.nsamples"])).astype(np.int64)
        gf = ak.to_numpy(ak.flatten(a["hycal.gain_factor"])).astype(np.float32)
        Ysrc = ak.to_numpy(ak.flatten(a["hycal.samples"], axis=1))
        N = mid_id.size
        st["n_events"] += ne
        st["n_chan"] += N
        midx = idx_of_id[mid_id]
        is_pw = np.zeros(N, dtype=bool)
        ok = midx >= 0
        is_pw[ok] = pw_arr[midx[ok]]
        if NS is None:
            u = nsm[is_pw]
            NS = int(np.bincount(u).argmax()) if u.size else int(wcfg.nsamples
                                                                 if hasattr(wcfg, "nsamples") else 100)
        mf_row = np.nonzero(is_pw & (nsm == NS))[0]
        P_h = np.zeros((N, MAX_PEAKS), np.float32)
        P_t = np.zeros((N, MAX_PEAKS), np.float32)
        P_i = np.zeros((N, MAX_PEAKS), np.float32)
        P_q = np.zeros((N, MAX_PEAKS), np.uint8)
        npk = np.zeros(N, np.uint8)
        kind = np.zeros(N, np.uint8)       # 1 mf, 2 wa_fb, 3 lg/other, 4 odd
        wa_mf = {}                         # row -> WA peaks (diag only)

        def put(r, plist, q):
            plist = plist[:MAX_PEAKS]
            npk[r] = len(plist)
            for k, (t, h, ii) in enumerate(plist):
                P_t[r, k], P_h[r, k], P_i[r, k], P_q[r, k] = t, h, ii, q

        if mf_row.size:
            res = rec.run(Ysrc[mf_row, :NS], midx[mf_row], overflow)
            st["n_sat_fallback"] += res.stats["n_sat_fallback"]
            pc = res.p_ch
            Eu = mf_pulse.pulse_energy_units(bank, midx[mf_row][pc], res.p_t0, res.p_amp, 1.0)
            keep = (res.p_amp > 0) & (Eu > 0)
            st["n_dropped"] += int((~keep).sum())
            pc, pt_, pa_, Eu = pc[keep], res.p_time[keep], res.p_amp[keep], Eu[keep]
            starts = np.searchsorted(pc, np.arange(mf_row.size + 1))
            kk = np.arange(pc.size) - starts[pc]
            fit = (kk < MAX_PEAKS) & ~res.wa_fb[pc]
            rr = mf_row[pc[fit]]
            P_t[rr, kk[fit]] = pt_[fit]
            P_h[rr, kk[fit]] = pa_[fit]
            P_i[rr, kk[fit]] = Eu[fit]
            P_q[rr, kk[fit]] = MF_BIT
            cnt = np.bincount(pc[fit], minlength=mf_row.size)
            npk[mf_row] = np.minimum(cnt, MAX_PEAKS)
            st["n_mf_pulses"] += int(fit.sum())
            fb = np.nonzero(res.wa_fb)[0]
            kind[mf_row] = 1
            for r in fb:
                row = mf_row[r]
                kind[row] = 2
                put(row, wa_peaks(wave_ana, Ysrc[row], NS), 0)
            if args.diag and start < args.diag_events:
                good = np.nonzero(fit & (pa_ >= 100.0))[0]
                rows_g = np.unique(mf_row[pc[good]])
                for row in rows_g:
                    wa_mf[row] = wa_peaks(wave_ana, Ysrc[row], NS)
                for g in good:
                    row = mf_row[pc[g]]
                    best = None
                    for (t, h, ii) in wa_mf[row]:
                        d = float(pt_[g]) - t
                        if abs(d) < 8.0 and (best is None or abs(d) < abs(best[0])) and ii > 0:
                            best = (d, ii)
                    if best is not None:
                        nm = names[midx[row]]
                        lst = ratio.setdefault(nm, ([], []))
                        lst[0].append(float(Eu[g]) / best[1])
                        lst[1].append(best[0])
        other = np.nonzero(kind == 0)[0]
        for row in other:
            if is_pw[row]:
                kind[row] = 4          # PbWO4 with unexpected sample count
            else:
                kind[row] = 3
            put(row, wa_peaks(wave_ana, Ysrc[row], int(nsm[row])), 0)
        has = npk > 0
        st["n_peak_chan"] += int(has.sum())
        st["n_mf"] += int(((kind == 1) & has).sum())
        st["n_wa_fb"] += int((kind == 2).sum())
        st["n_wa_fb_peaks"] += int(((kind == 2) & has).sum())
        st["n_odd"] += int((kind == 4).sum())
        st["n_odd_peaks"] += int(((kind == 4) & has).sum())
        st["n_lg"] += int((kind == 3).sum())
        st["n_lg_peaks"] += int(((kind == 3) & has).sum())

        off = np.concatenate([[0], np.cumsum(nch)])
        off8 = off * MAX_PEAKS

        def jag(v, o):
            return ak.unflatten(v, np.diff(o))

        flat["flat"].extend(dict(
            event_num=evn, trigger_type=ttype, trigger_bits=tbits, nch=nch,
            n8=(nch * MAX_PEAKS).astype(np.int32),
            module_id=jag(mid_id.astype(np.uint16), off), module_type=jag(mtype, off),
            gain_factor=jag(gf, off), npeaks=jag(npk, off),
            pk_h=jag(P_h.reshape(-1), off8), pk_t=jag(P_t.reshape(-1), off8),
            pk_i=jag(P_i.reshape(-1), off8), pk_q=jag(P_q.reshape(-1), off8)))
        print(f"[writer] {stop}/{n_entries} entries, {st['n_events']} events written "
              f"({time.time() - t_start:.0f} s)", flush=True)
    flat.close()

    macro = os.path.join(tmpd, "fill_peaks.C")
    with open(macro, "w") as fh:
        fh.write(MACRO)
    tmp_out = os.path.join(tmpd, "out.root")
    r = subprocess.run(["root", "-l", "-b", "-q", f'fill_peaks.C+("{flat_path}","{tmp_out}")'],
                       cwd=tmpd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(tmp_out):
        sys.stderr.write(r.stdout + r.stderr)
        sys.exit("[ERROR] root macro failed (tmp dir kept: %s)" % tmpd)
    shutil.move(tmp_out, out_path)
    shutil.rmtree(tmpd, ignore_errors=True)

    meta = dict(st, in_path=in_path, in_sha1=sha1_file(in_path), bank=args.mf_templates,
                bank_sha1=info["sha1"], bank_fallback=len(info["fallback"]),
                argv=sys.argv, out=out_path, kappa=1.0, ns=NS,
                seconds=round(time.time() - t_start, 1))
    with open(os.path.join(args.out, stem + "_writer.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    if args.diag:
        with open(os.path.join(args.out, stem + "_ratio.tsv"), "w") as fh:
            fh.write("module\tn_match\tr_med\tr_iqr\tdt_med_ns\n")
            for nm in sorted(ratio):
                r_, d_ = (np.array(x) for x in ratio[nm])
                q1, q3 = np.percentile(r_, [25, 75])
                fh.write(f"{nm}\t{r_.size}\t{np.median(r_):.5f}\t{q3 - q1:.5f}\t"
                         f"{np.median(d_):.4f}\n")
    print(json.dumps({k: v for k, v in meta.items() if k.startswith("n_")}), flush=True)


if __name__ == "__main__":
    main()
