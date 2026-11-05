import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 9. Which chromosome arms are gained or lost?

    **Lecture 3.** For each tumour and chromosome arm: the median log2 copy
    ratio of the arm's genes. An arm counts as **gained** above +0.2 and
    **lost** below −0.2. (The p arms of chromosomes 13, 14, 15, 21 and 22
    carry almost no genes and are left out.)

    Some changes recur in particular cancers: +7 and −10 in glioblastoma,
    −3p in kidney cancer, +8q in many carcinomas. Are they chosen, or just
    tolerated?
    """)
    return


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    import tcga

    tcga.style()
    return mo, np, plt, tcga


@app.cell
def _(mo, tcga):
    ARM_LOG2 = 0.2
    arms = tcga.load("arms").merge(tcga.load("patients"))
    _any = arms.groupby("patient_id").log2_ratio.apply(
        lambda v: (v.abs() > ARM_LOG2).any()
    )
    mo.md(
        f"**{100 * _any.mean():.0f}%** of {len(_any):,} tumours have at least "
        "one arm gained or lost."
    )
    return ARM_LOG2, arms


@app.cell
def _(ARM_LOG2, arms, tcga):
    TYPES = [
        "OV", "ESCA", "LUSC", "HNSC", "PAAD", "COADREAD", "STAD", "BLCA",
        "LUAD", "BRCA", "UCEC", "LIHC", "GBM", "SARC", "SKCM", "PRAD", "KIRC",
        "THCA",
    ]  # fmt: skip
    ORDER = [
        f"{n}{a}" for n in range(1, 23) for a in "pq"
        if not (n in {13, 14, 15, 21, 22} and a == "p")
    ]  # fmt: skip
    _a = arms[arms.cancer_type.isin(TYPES)].assign(
        gain=lambda d: d.log2_ratio > ARM_LOG2,
        loss=lambda d: d.log2_ratio < -ARM_LOG2,
    )
    _g = _a.groupby(["cancer_type", "arm"])
    gain = _g.gain.mean().unstack().loc[TYPES, ORDER]
    loss = _g.loss.mean().unstack().loc[TYPES, ORDER]
    names = [tcga.CANCER_NAMES[t] for t in TYPES]
    return ORDER, gain, loss, names


@app.cell
def _(ORDER, gain, loss, names, np, plt):
    # gains up (red), losses down (blue): one diverging scale
    _net = gain.to_numpy() - loss.to_numpy()
    fig, ax = plt.subplots(figsize=(11, 5))
    _im = ax.imshow(100 * _net, cmap="RdBu_r", vmin=-60, vmax=60,
                    aspect="auto")  # fmt: skip
    ax.set_xticks(np.arange(len(ORDER)), ORDER, rotation=90)
    ax.set_yticks(np.arange(len(names)), names)
    ax.spines[:].set_visible(False)
    fig.colorbar(_im, ax=ax, shrink=0.7,
                 label="% gained minus % lost")  # fmt: skip
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. In glioblastoma, how many tumours have both +7 and −10? Compare with
       the number expected if the two were independent.
    2. 8q carries MYC. Is 8q gained more often than 8p? What does 8p loss
       plus 8q gain suggest about how the change arose (an isochromosome)?
    """)
    return


if __name__ == "__main__":
    app.run()
