#!/usr/bin/env python3
"""
Post-processing of a Turek-Hron run

Input: displacement of point A (lines "TipA t ux uy" in the solver log) and the forces
acting on the cylinder and the flag (postProcessing/forcesTurek/*/force.dat).
For the final time window the mean value, amplitude and frequency of the oscillations are
evaluated as in the benchmark: mean = (max + min)/2, amplitude = (max - min)/2; the
frequency is obtained from upward crossings of the mean level. The results are compared
with the reference values and plotted to turek_results.png

    python3 scripts/postprocess.py                              # averaging window of 2 s
    python3 scripts/postprocess.py --window 3 --log log.pimpleFoam
"""
import argparse
import re
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent.parent

# Turek & Hron 2006: (mean, amplitude, frequency); displacements in m, forces in N/m
REF = {
    "FSI2": dict(ux=(-14.58e-3, 12.44e-3, 3.86), uy=(1.23e-3, 80.6e-3, 1.93),
                 drag=(208.83, 73.75, 3.86), lift=(0.88, 234.2, 1.93)),
    "FSI3": dict(ux=(-2.69e-3, 2.53e-3, 10.9), uy=(1.48e-3, 34.38e-3, 5.46),
                 drag=(457.3, 22.66, 10.9), lift=(2.22, 149.78, 5.46)),
}


def read_tip(log):
    t, ux, uy = [], [], []
    for line in Path(log).read_text(errors="ignore").splitlines():
        if line.startswith("TipA "):
            a = line.split()
            t.append(float(a[1])); ux.append(float(a[2])); uy.append(float(a[3]))
    return np.array(t), np.array(ux), np.array(uy)


def read_forces():
    files = sorted((CASE/"postProcessing"/"forcesTurek").glob("*/force.dat"),
                   key=lambda p: float(p.parent.name))
    rows = []
    for f in files:
        for line in f.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            nums = [float(x) for x in re.sub(r"[()]", " ", line).split()]
            if len(nums) >= 10:
                rows.append((nums[0], nums[1], nums[2]))
            elif len(nums) >= 7:
                rows.append((nums[0], nums[1] + nums[4], nums[2] + nums[5]))
    if not rows:
        return None
    a = np.array(rows)
    _, idx = np.unique(a[:, 0], return_index=True)
    return a[idx].T


def stats(t, u, window):
    m = t >= t[-1] - window
    t, u = t[m], u[m]
    mean, amp = 0.5*(u.max() + u.min()), 0.5*(u.max() - u.min())
    s = np.sign(u - mean)
    up = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    freq = (len(up) - 1)/(t[up[-1]] - t[up[0]]) if len(up) > 2 else float("nan")
    return mean, amp, freq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", default=str(CASE/"log.pimpleFoam"))
    ap.add_argument("--window", type=float, default=2.0)
    ap.add_argument("--case", default=None)
    a = ap.parse_args()

    case = a.case
    if case is None and (CASE/"caseInfo.txt").exists():
        case = (CASE/"caseInfo.txt").read_text().split()[1]
    ref = REF.get(case, {})

    series = {}
    t, ux, uy = read_tip(a.log)
    if len(t):
        series["ux"], series["uy"] = (t, ux), (t, uy)
    F = read_forces()
    if F is not None:
        series["drag"], series["lift"] = (F[0], F[1]), (F[0], F[2])

    print(f"Case {case}, averaging window: last {a.window} s")
    print(f"{'':6s}{'computed: mean ± amp [Hz]':>36s}   {'reference: mean ± amp [Hz]':>36s}")
    for key in ("ux", "uy", "drag", "lift"):
        if key not in series:
            continue
        mean, amp, f = stats(*series[key], a.window)
        sc, unit = (1e3, "e-3 m") if key in ("ux", "uy") else (1.0, "N/m  ")
        line = f"{key:6s}{mean*sc:12.3f} ± {amp*sc:9.3f} [{f:6.3f}] {unit}"
        if key in ref:
            r = ref[key]
            line += f"   {r[0]*sc:12.3f} ± {r[1]*sc:9.3f} [{r[2]:6.3f}]"
        print(line)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = [k for k in ("ux", "uy", "drag", "lift") if k in series]
        fig, axs = plt.subplots(len(keys), 1, figsize=(11, 2.4*len(keys)), sharex=True)
        for ax, k in zip(np.atleast_1d(axs), keys):
            tt, u = series[k]
            sc = 1e3 if k in ("ux", "uy") else 1.0
            ax.plot(tt, u*sc, lw=0.8)
            if k in ref and ref[k][1] > 0:
                for v in (ref[k][0] - ref[k][1], ref[k][0] + ref[k][1]):
                    ax.axhline(v*sc, color="tab:red", ls="--", lw=0.8)
            ax.set_ylabel(k + (" ×1e-3 m" if sc > 1 else ", N/m"))
            ax.grid(alpha=0.3)
        np.atleast_1d(axs)[-1].set_xlabel("t, s")
        fig.suptitle(f"Turek-Hron {case}: computed (solid), reference extrema (dashed)")
        fig.tight_layout()
        fig.savefig(CASE/"turek_results.png", dpi=150)
        print("Plots saved to turek_results.png")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
