import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 10. Extra copies, extra RNA: gene dosage on 8q

    **Lecture 3.** Does gaining a chromosome arm change the expression of
    every gene on it, or only a few drivers?

    Breast tumours with 8q **gained** (median log2 ratio above +0.2) against
    tumours with 8q **neutral** (between −0.1 and +0.1). For each
    expressed gene (mean log2 TPM above 1): the mean expression difference
    between the two groups, for genes on 8q against 1,500 genes elsewhere.

    If most 8q genes go up, an arm gain raises hundreds of proteins at
    once: the cell has to fold, assemble or degrade the extra protein
    (proteotoxic stress).
    """)
    return


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    import tcga
    from colors import COLOR

    tcga.style()
    return COLOR, mo, np, plt, tcga


@app.cell
def _(tcga):
    _arms = tcga.load("arms")
    _8q = _arms[_arms.arm == "8q"].set_index("patient_id").log2_ratio
    group = _8q.map(
        lambda v: "gain" if v > 0.2 else "neutral" if abs(v) < 0.1 else None
    ).dropna()

    expr = tcga.load("expression_brca")
    expr = expr[expr.patient_id.isin(group.index)]
    expr["group"] = expr.patient_id.map(group)
    expr.group.value_counts()
    return (expr,)


@app.cell
def _(expr, tcga):
    _genes = tcga.load("genes").set_index("gene_id")
    _cols = [c for c in expr.columns if c.startswith("ENSG")]
    _means = expr.groupby("group")[_cols].mean().T
    change = _means.assign(
        diff=_means["gain"] - _means["neutral"],
        mean=expr[_cols].mean(),
        on8q=_genes.arm.reindex(_cols).eq("8q").to_numpy(),
        symbol=_genes.gene_symbol.reindex(_cols).to_numpy(),
    )
    change = change[change["mean"] > 1]  # expressed genes only
    change.sort_values("diff", ascending=False).head(10)
    return (change,)


@app.cell
def _(COLOR, change, mo, np, plt):
    _bins = np.arange(-1.5, 2.01, 0.1)
    _a, _b = change[change.on8q], change[~change.on8q]
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.hist(_b["diff"], bins=_bins, density=True, color=COLOR.GREY.value,
            label=f"other genes (n={len(_b)})")  # fmt: skip
    ax.hist(_a["diff"], bins=_bins, density=True, color=COLOR.PINK.value,
            alpha=0.8, label=f"8q genes (n={len(_a)})")  # fmt: skip
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("expression change, 8q gained vs neutral (log2)")
    ax.set_ylabel("density")
    ax.legend(loc="upper right")
    fig.tight_layout()
    mo.vstack([
        fig,
        mo.md(
            f"8q genes: median change **{_a['diff'].median():+.2f}** log2, "
            f"**{100 * (_a['diff'] > 0).mean():.0f}%** up. Other genes: "
            f"{_b['diff'].median():+.2f}."
        ),
    ])  # fmt: skip
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Is MYC (ENSG00000136997) among the most increased 8q genes? Look it
       up in `change`.
    2. One copy more on a diploid background is 3/2 copies, log2(1.5) =
       0.58. How close is the median 8q gene to that? What could make it
       smaller (tumour purity, partial gains, dosage compensation)?
    """)
    return


if __name__ == "__main__":
    app.run()
