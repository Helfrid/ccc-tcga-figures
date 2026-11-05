import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 4. Breaking the Rb pathway

    **Lecture 2.** p16 (CDKN2A) inhibits cyclin D–CDK4/6, which inactivates
    Rb. A tumour can release the brake at any step:

    - **lose** a brake: CDKN2A or RB1 deleted on both copies, or mutated;
    - **gain** an accelerator: CCND1, CDK4 or CDK6 amplified.

    Most tumours use one route. Losing p16 *and* Rb together is pointless:
    without Rb, p16 has nothing left to protect.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    import tcga
    from colors import COLOR

    tcga.style()
    return COLOR, mo, np, pd, plt, tcga


@app.cell
def _(tcga):
    GENES = ["CDKN2A", "CCND1", "CDK4", "CDK6", "RB1"]
    SUPPRESSORS = {"CDKN2A", "RB1"}

    mutations = tcga.load("mutations")
    cn = tcga.copy_number_wide()
    _p = tcga.sequenced()
    _p = _p[_p.patient_id.isin(cn.index)]
    ctype = dict(zip(_p.patient_id, _p.cancer_type))

    # per gene: route -> set of patients
    routes = {}
    for _g in GENES:
        if _g in SUPPRESSORS:
            routes[_g] = {
                "del": tcga.deleted(cn, _g),
                "mut": tcga.inactivated(mutations, _g),
            }
        else:
            routes[_g] = {"amp": tcga.amplified(cn, _g)}
    altered = {
        g: set().union(*r.values()) & set(ctype) for g, r in routes.items()
    }
    return GENES, altered, ctype, routes


@app.cell
def _(GENES, altered, ctype, pd, tcga):
    TYPES = [
        "GBM", "SKCM", "BLCA", "ESCA", "HNSC", "LUSC",
        "PAAD", "SARC", "BRCA", "LUAD", "STAD", "COADREAD",
    ]  # fmt: skip
    _rows = {}
    for _t in TYPES:
        _pts = {p for p, t in ctype.items() if t == _t}
        _any = set().union(*(altered[g] & _pts for g in GENES))
        _rows[tcga.CANCER_NAMES[_t]] = {
            g: 100 * len(altered[g] & _pts) / len(_pts) for g in GENES
        } | {"any": 100 * len(_any) / len(_pts), "n": len(_pts)}
    landscape = pd.DataFrame(_rows).T
    landscape.round(1)
    return (landscape,)


@app.cell
def _(GENES, landscape, np, plt):
    _cols = GENES + ["any"]
    _v = landscape[_cols].to_numpy()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.imshow(_v, cmap="Blues", vmin=0, vmax=80, aspect="auto")
    ax.set_xticks(range(len(_cols)), _cols)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(landscape)), landscape.index)
    for (_i, _j), _x in np.ndenumerate(_v):
        if _x >= 1:
            ax.text(_j, _i, f"{_x:.0f}", ha="center", va="center", fontsize=7,
                    color="white" if _x > 45 else "black")  # fmt: skip
    ax.axvline(len(GENES) - 0.5, color="white", lw=3)
    ax.spines[:].set_visible(False)
    ax.set_title("% of tumours with each gene altered", loc="left", pad=24)
    fig.tight_layout()
    fig
    return


@app.cell
def _(altered, ctype, mo):
    # CDKN2A and RB1 together, across all cancers
    _pan = set(ctype)
    _a = altered["CDKN2A"] & _pan
    _b = altered["RB1"] & _pan
    _expected = len(_a) * len(_b) / len(_pan)
    mo.md(
        f"Across {len(_pan):,} tumours: CDKN2A lost in {len(_a):,}, RB1 lost "
        f"in {len(_b):,}. Both: **{len(_a & _b)}** observed against "
        f"**{_expected:.1f}** expected if independent."
    )
    return


@app.cell
def _(COLOR, GENES, ctype, np, plt, routes):
    # one tumour per column, glioblastoma
    _pts = {p for p, t in ctype.items() if t == "GBM"}
    _calls = {p: {} for p in _pts}
    for _g in GENES:
        for _kind, _set in routes[_g].items():
            for _pt in _set & _pts:
                _calls[_pt][_g] = "both" if _g in _calls[_pt] else _kind
    _hit = sorted(
        (p for p in _pts if _calls[p]),
        key=lambda p: tuple(0 if g in _calls[p] else 1 for g in GENES),
    )
    _colour = {"del": COLOR.BLUE.value, "mut": COLOR.YELLOW.value,
               "amp": COLOR.PINK.value, "both": COLOR.PURPLE.value}  # fmt: skip
    fig2, ax2 = plt.subplots(figsize=(9, 2.9))
    for _j, _p in enumerate(_hit):
        for _i, _g in enumerate(GENES):
            _c = _calls[_p].get(_g)
            ax2.add_patch(plt.Rectangle(
                (_j, _i - 0.4), 0.8, 0.8, lw=0,
                color=_colour[_c] if _c else COLOR.GREY.value,
                alpha=1 if _c else 0.35,
            ))  # fmt: skip
    ax2.set_xlim(-1, len(_hit) + 1)
    ax2.set_ylim(len(GENES) - 0.5, -0.5)
    ax2.set_yticks(np.arange(len(GENES)), GENES)
    ax2.set_xticks([])
    ax2.spines[:].set_visible(False)
    ax2.set_xlabel(f"Glioblastoma: {len(_hit)} of {len(_pts)} tumours altered")
    for _k, _lab in [("del", "deleted"), ("mut", "mutated"),
                     ("amp", "amplified"), ("both", "both")]:  # fmt: skip
        ax2.scatter([], [], marker="s", color=_colour[_k], label=_lab)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4)
    fig2.tight_layout()
    fig2
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Which cancer type relies most on RB1 loss rather than CDKN2A loss?
       Is CDKN2A lost in the same cancer types?
    2. CDKN2A sits next to CDKN2B and MTAP on chromosome 9p21. Look at
       `cn["CDKN2B"]`: how often is it deleted together with CDKN2A?
    """)
    return


if __name__ == "__main__":
    app.run()
