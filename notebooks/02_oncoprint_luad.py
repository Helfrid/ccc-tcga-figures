import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 2. One driver per tumour: lung adenocarcinoma

    **Lecture 2.** An *oncoprint*: one column per tumour, one row per gene,
    a mark where the gene is altered. In lung adenocarcinoma KRAS, EGFR,
    BRAF and ERBB2 are each common, yet a tumour almost never has two of
    them. One hit to the same pathway is enough; a second one adds nothing
    (and too much RAS signalling can even stop cells, oncogene-induced
    senescence).

    We compare the number of tumours with two or more hits against the
    number expected if the genes were hit independently.
    """)
    return


@app.cell
def _():
    import math

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    import tcga
    from colors import COLOR

    tcga.style()
    return COLOR, math, mo, np, plt, tcga


@app.cell
def _(tcga):
    GENES = ["KRAS", "EGFR", "BRAF", "ERBB2", "NRAS"]
    CANCER = "LUAD"

    mutations = tcga.load("mutations")
    cn = tcga.copy_number_wide()
    _p = tcga.sequenced()
    pts = set(_p.patient_id[(_p.cancer_type == CANCER)]) & set(cn.index)

    # per patient and gene: "mut", "amp" or "both"
    calls = {p: {} for p in pts}
    for _g in GENES:
        _mut = tcga.mutated(mutations, _g) & pts
        _amp = (
            tcga.amplified(cn, _g) & pts if _g in ("EGFR", "ERBB2") else set()
        )
        for _pt in _mut | _amp:
            calls[_pt][_g] = (
                "both" if _pt in _mut and _pt in _amp
                else "mut" if _pt in _mut else "amp"
            )  # fmt: skip
    len(pts)
    return GENES, calls, pts


@app.cell
def _(GENES, calls, math, mo, pts):
    n = len(pts)
    freq = {g: sum(g in calls[p] for p in pts) / n for g in GENES}
    observed = sum(len(c) >= 2 for c in calls.values())

    # P(no gene hit) and P(exactly one hit) if genes were independent
    _none = math.prod(1 - freq[g] for g in GENES)
    _one = sum(
        freq[g] * math.prod(1 - freq[h] for h in GENES if h != g)
        for g in GENES
    )
    expected = n * (1 - _none - _one)
    mo.md(
        f"**{n} tumours.** "
        + ", ".join(f"{g} {100 * f:.1f}%" for g, f in freq.items())
        + f". Two or more hits: **{observed}** observed, "
        f"**{expected:.1f}** expected if independent."
    )
    return (freq,)


@app.cell
def _(COLOR, GENES, calls, freq, np, plt):
    # altered tumours only, sorted gene by gene (the classic staircase)
    hit = [p for p in calls if calls[p]]
    hit.sort(key=lambda p: tuple(0 if g in calls[p] else 1 for g in GENES))
    _colour = {"mut": COLOR.BLUE.value, "amp": COLOR.PINK.value,
               "both": COLOR.PURPLE.value}  # fmt: skip

    fig, ax = plt.subplots(figsize=(9, 3))
    for _j, _p in enumerate(hit):
        for _i, _g in enumerate(GENES):
            _c = calls[_p].get(_g)
            ax.add_patch(
                plt.Rectangle(
                    (_j, _i - 0.4), 0.8, 0.8,
                    color=_colour[_c] if _c else COLOR.GREY.value,
                    alpha=1 if _c else 0.35, lw=0,
                )
            )  # fmt: skip
    ax.set_xlim(-1, len(hit) + 1)
    ax.set_ylim(len(GENES) - 0.5, -0.5)
    ax.set_yticks(
        np.arange(len(GENES)),
        [f"{g}  {100 * freq[g]:.0f}%" for g in GENES],
    )
    ax.set_xticks([])
    ax.spines[:].set_visible(False)
    ax.set_xlabel(f"{len(hit)} altered tumours, one column each")
    for _k, _lab in [("mut", "mutated"), ("amp", "amplified"),
                     ("both", "both")]:  # fmt: skip
        ax.scatter([], [], marker="s", color=_colour[_k], label=_lab)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3)
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Change `CANCER` to `"COADREAD"` (colorectal). Are KRAS and BRAF still
       exclusive? What about KRAS and NRAS?
    2. Add `"PIK3CA"` to `GENES`, after defining its hotspots with
       `tcga.HOTSPOTS["PIK3CA"] = {542, 545, 1047}`. Does it avoid KRAS,
       or co-occur? Why might that be?
    """)
    return


if __name__ == "__main__":
    app.run()
