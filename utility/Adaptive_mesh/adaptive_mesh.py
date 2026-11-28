#!/usr/bin/env python3
"""
adaptive_mesh.py  -  generate an x-grid that is dense near peaks / steep regions
                     of a curve y(x) and sparse where the curve is flat.

USAGE
    python adaptive_mesh.py  datafile  [key=value ...]

EXAMPLES
    python adaptive_mesh.py energy.dat x=1 y=3
    python adaptive_mesh.py energy.dat x=1 y=3 err=5 N=30
    python adaptive_mesh.py energy.dat x=1 y=3 N=20 A=50 p=3 kc=2
    python adaptive_mesh.py energy.dat x = 1 y = 3 N = 25     (spaces allowed)

PARAMETERS  (all optional; defaults are tuned to put many points at the peak)
    x    column number of x (1-based)                                  [1]
    y    column number of y (1-based)                                  [2]
    err  column number of y-errors (1-based); 0 = none, noise estimated [0]
    N    total number of grid points                                   [30]
    A    contrast: peak density up to (1+A) x flat-region density       [30]
         larger -> more points at the peak;  A=0 -> uniform grid
    p    sharpness: >1 squeezes the dense region onto the peak         [2]
         larger -> narrower, denser cluster;  p<1 -> broader cluster
    kc   curvature weight vs slope                                     [1]
         larger -> points on the peak TOP;  0 -> only on the steep sides
    s    smoothing factor of the spline (x number of data points)      [1]
         smaller -> follows data more closely (and its noise)
    dxmin  minimum allowed spacing between points (0 = no limit)       [0]
    dxmax  maximum allowed spacing between points (0 = no limit)       [0]
    xmin, xmax   range of the new grid      [data range]
    dec  decimals for the output x values   [5]
    order  order of the points in the output: desc (high -> low) or asc [desc]
    out  output file name (without extension) [grid]

OUTPUT
    <out>.out   line 1: "# number of points", line 2: "# all points in one line",
                then two columns: new x points and smoothed y at those points
    <out>.png   plot: data, smoothed curve, new points, spacing
"""
import re
import sys

import numpy as np
from scipy.interpolate import UnivariateSpline

_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy 1.x / 2.x

DEFAULTS = dict(x=1, y=2, err=0, N=30, A=30.0, p=2.0, kc=1.0, s=1.0,
                dxmin=0.0, dxmax=0.0,
                xmin=None, xmax=None, dec=5, out="grid", order="desc")
INT_KEYS = {"x", "y", "err", "N", "dec"}


def parse_args(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0)
    fname = argv[1]
    text = " ".join(argv[2:])
    pairs = re.findall(r"(\w+)\s*=\s*([^\s=]+)", text)
    leftover = re.sub(r"(\w+)\s*=\s*([^\s=]+)", "", text).strip()
    if leftover:
        sys.exit(f"Could not understand: '{leftover}'. Use key=value, e.g. N=30")
    par = dict(DEFAULTS)
    lower = {k.lower(): k for k in DEFAULTS}
    for k, v in pairs:
        key = k if k in DEFAULTS else lower.get(k.lower())
        if key is None:
            sys.exit(f"Unknown parameter '{k}'. Allowed: {', '.join(DEFAULTS)}")
        if key == "order":
            v = v.lower()
            if v not in ("desc", "asc"):
                sys.exit("order must be 'desc' or 'asc'")
            par[key] = v
        elif key == "out":
            par[key] = v
        elif key in INT_KEYS:
            par[key] = int(v)
        else:
            par[key] = float(v)
    return fname, par


def build_grid(X, Y, E, par):
    # sort and merge duplicate x values
    o = np.argsort(X)
    X, Y = X[o], Y[o]
    E = E[o] if E is not None else None
    X, idx = np.unique(X, return_index=True)
    Y = Y[idx]
    E = E[idx] if E is not None else None
    n = len(X)
    if n < 6:
        sys.exit("Need at least 6 data points.")

    # weights for the smoothing spline
    if E is None or np.any(E <= 0):
        # noise estimate from second differences (rms -> conservative, so
        # the spline does not chase noise near the peak, where MC noise is largest)
        d2 = np.diff(Y, 2)
        sigma = np.sqrt(np.mean(d2**2) / 6.0)
        if sigma <= 0:
            sigma = 1e-6 * (np.ptp(Y) or 1.0)
        E = np.full(n, sigma)
    spl = UnivariateSpline(X, Y, w=1.0 / E, k=4, s=par["s"] * n)

    a = X.min() if par["xmin"] is None else par["xmin"]
    b = X.max() if par["xmax"] is None else par["xmax"]
    xf = np.linspace(a, b, 20001)
    c = spl(xf)
    scale = np.ptp(c) or 1.0
    L = b - a
    d1 = spl.derivative(1)(xf) * L / scale            # dimensionless slope
    d2 = spl.derivative(2)(xf) * L**2 / scale         # dimensionless curvature

    g = np.sqrt(d1**2 + par["kc"] * np.abs(d2))       # feature strength
    w = 1.0 + par["A"] * (g / (g.max() or 1.0)) ** par["p"]   # point density

    # Spacing limits: local spacing ~ mean(w)*L/((N-1)*w), so dx >= dxmin
    # <=> w <= mean(w)*L/((N-1)*dxmin). Clip w and iterate (mean(w) changes).
    M = par["N"] - 1
    if par["dxmin"] > 0 and par["dxmin"] * M > L:
        sys.exit(f"dxmin={par['dxmin']} too large: N={par['N']} points need "
                 f"dxmin <= {L / M:.4g}")
    if par["dxmax"] > 0 and par["dxmax"] * M < L:
        sys.exit(f"dxmax={par['dxmax']} too small: N={par['N']} points need "
                 f"dxmax >= {L / M:.4g}")
    w0 = w.copy()
    for _ in range(200):
        wm = _trapz(w, xf) / L
        hi = wm * L / (M * par["dxmin"]) if par["dxmin"] > 0 else np.inf
        lo = wm * L / (M * par["dxmax"]) if par["dxmax"] > 0 else 0.0
        w_new = np.clip(w0, lo, hi)
        if np.allclose(w_new, w, rtol=1e-10):
            break
        w = w_new

    W = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(xf))])
    x_new = np.interp(np.linspace(0, 1, par["N"]), W / W[-1], xf)
    x_new = np.unique(np.round(x_new, par["dec"]))
    return X, Y, E, spl, xf, x_new


def make_plot(X, Y, E, spl, xf, x_new, par, xlab, ylab, has_err):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # compact figure with large fonts (small file, readable in a README)
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 14,
                         "legend.fontsize": 11.5, "xtick.labelsize": 12,
                         "ytick.labelsize": 12})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1.2]})
    if has_err:
        ax1.errorbar(X, Y, E, fmt="o", ms=3.5, mfc="none", color="0.6",
                     capsize=1.5, label="input data")
    else:
        ax1.plot(X, Y, "o", ms=3.5, mfc="none", color="0.6", label="input data")
    ax1.plot(xf, spl(xf), color="#c0504d", lw=1.8, label="smoothed curve",
             zorder=2)
    ax1.plot(x_new, spl(x_new), "o", ms=6, color="#1f4e79", zorder=3,
             label=f"adaptive grid ({len(x_new)} pts)")
    lo = ax1.get_ylim()[0]
    ax1.plot(x_new, np.full_like(x_new, lo), "|", ms=12, color="#1f4e79")
    ax1.set_ylabel(ylab)
    ax1.legend(frameon=False, loc="best")
    ax1.set_title(f"N={par['N']}  A={par['A']:g}  p={par['p']:g}  "
                  f"kc={par['kc']:g}  s={par['s']:g}  dxmin={par['dxmin']:g}  "
                  f"dxmax={par['dxmax']:g}", fontsize=10.5)
    ax2.plot(0.5 * (x_new[1:] + x_new[:-1]), np.diff(x_new), "o-", ms=3,
             color="#1f4e79")
    ax2.set_ylabel(r"$\Delta x$")
    ax2.set_xlabel(xlab)
    fig.tight_layout()
    png = f"{par['out']}.png"
    fig.savefig(png, dpi=120, bbox_inches="tight", pad_inches=0.05)
    try:  # shrink the file: 128-colour palette PNG (looks the same)
        from PIL import Image
        Image.open(png).convert("RGB").quantize(
            colors=128, method=Image.Quantize.MEDIANCUT).save(png, optimize=True)
    except Exception:
        pass


def main():
    fname, par = parse_args(sys.argv)
    data = np.loadtxt(fname, comments="#", ndmin=2)
    ncol = data.shape[1]
    for k in ("x", "y") + (("err",) if par["err"] else ()):
        if not 1 <= par[k] <= ncol:
            sys.exit(f"{k}={par[k]} is out of range: the file has {ncol} columns.")

    X = data[:, par["x"] - 1]
    Y = data[:, par["y"] - 1]
    E = data[:, par["err"] - 1] if par["err"] else None

    X, Y, E, spl, xf, x_new = build_grid(X, Y, E, par)

    x_w = x_new[::-1] if par["order"] == "desc" else x_new   # output order
    out = f"{par['out']}.out"
    # line 1: number of points; line 2: all points on one line; then x, y rows
    hdr = (f"{len(x_new)}\n"
           + " ".join(f"{v:.{par['dec']}f}" for v in x_w))
    np.savetxt(out, np.column_stack([x_w, spl(x_w)]),
               fmt=[f"%12.{par['dec']}f", "%14.6e"], header=hdr)
    make_plot(X, Y, E, spl, xf, x_new, par, f"column {par['x']}",
              f"column {par['y']}", par["err"] > 0)

    dx = np.diff(x_new)
    ipk = np.argmax(spl(xf))
    print(f"{len(x_new)} points written to {out} (plot: {par['out']}.png)")
    print(f"smoothed peak at x = {xf[ipk]:.{par['dec']}f}")
    print(f"spacing: min {dx.min():.3g} at x~{x_new[np.argmin(dx)]:.{par['dec']}f}, "
          f"max {dx.max():.3g}")
    print(" ".join(f"{v:.{par['dec']}f}" for v in x_w))


if __name__ == "__main__":
    main()
