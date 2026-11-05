import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 7. Without p53, genomes come apart

    **Lecture 3.** For each tumour, the percent of genes with a changed copy
    number: about 3,600 genes along the autosomes, counted as changed when
    the log2 copy ratio is beyond ±0.3 (roughly a quarter more or a fifth
    less DNA than the tumour's own average).

    Because it is measured against the tumour's own average, a whole-genome
    doubling alone does not count (that is notebook 8).

    Each row: one cancer type with at least 20 TP53-mutant and 20 wild-type
    tumours. Dots: the median tumour in each group.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt

    import tcga
    from colors import COLOR

    tcga.style()
    return COLOR, mo, plt, tcga


@app.cell
def _(tcga):
    _m = tcga.load("mutations")
    _coding = tcga.MISSENSE + tcga.TRUNCATING
    _tp53 = set(
        _m.patient_id[
            (_m.gene_symbol == "TP53")
            & _m.variant_classification.isin(_coding)
        ]
    )
    tumours = tcga.sequenced().merge(tcga.load("genome_altered"))
    tumours["tp53"] = tumours.patient_id.isin(_tp53)
    tumours.head()
    return (tumours,)


@app.cell
def _(mo, tcga, tumours):
    _g = tumours.groupby(["cancer_type", "tp53"]).fraction_altered
    _t = _g.median().unstack().mul(100)
    _n = _g.size().unstack()
    by_type = (
        _t[(_n[True] >= 20) & (_n[False] >= 20)]
        .rename(columns={True: "mutant", False: "wild type"})
        .sort_values("mutant")
        .rename(index=tcga.CANCER_NAMES)
    )
    _higher = (by_type["mutant"] > by_type["wild type"]).sum()
    _all = tumours.groupby("tp53").fraction_altered.median().mul(100)
    mo.md(
        f"Higher in TP53-mutant tumours in **{_higher} of {len(by_type)}** "
        f"cancer types. All tumours: {_all[False]:.0f}% of genes changed "
        f"(wild type) against {_all[True]:.0f}% (mutant)."
    )
    return (by_type,)


@app.cell
def _(COLOR, by_type, plt):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    _y = range(len(by_type))
    ax.hlines(_y, by_type["wild type"], by_type["mutant"],
              color=COLOR.GREY.value, lw=3)  # fmt: skip
    ax.scatter(by_type["wild type"], _y, color=COLOR.LIGHT_BLUE.value, s=40,
               label="TP53 wild type", zorder=3)  # fmt: skip
    ax.scatter(by_type["mutant"], _y, color=COLOR.PINK.value, s=40,
               label="TP53 mutant", zorder=3)  # fmt: skip
    ax.set_yticks(_y, by_type.index)
    ax.set_xlabel("% of genes with a changed copy number (median tumour)")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Do point mutations follow the same pattern? Count mutations per
       patient and compare TP53-mutant and wild-type tumours within each
       type. (Look at endometrium and colorectal: why are they different?)
    2. Plot the whole distribution for breast cancer rather than the
       median (a histogram or a box plot per group).
    """)
    return


if __name__ == "__main__":
    app.run()
