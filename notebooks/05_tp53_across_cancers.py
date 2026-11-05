import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 5. TP53: the gene cancer loses most often

    **Lecture 3.** The percent of tumours with a TP53 mutation (any
    missense, in-frame or truncating change), in each cancer type.
    Ovarian cancer is left out (too few mutation calls in TCGA); high-grade
    serous ovarian cancer is close to 100% TP53 mutant.
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
    patients = tcga.sequenced()
    _m = tcga.load("mutations")
    _coding = tcga.MISSENSE + tcga.TRUNCATING
    tp53 = set(
        _m.patient_id[
            (_m.gene_symbol == "TP53")
            & _m.variant_classification.isin(_coding)
        ]
    )
    patients = patients.assign(tp53=patients.patient_id.isin(tp53))
    return (patients,)


@app.cell
def _(mo, patients, tcga):
    TYPES = [
        "ESCA", "LUSC", "HNSC", "PAAD", "COADREAD", "STAD", "BLCA", "LUAD",
        "BRCA", "UCEC", "LIHC", "GBM", "SARC", "SKCM", "PRAD", "KIRC", "THCA",
    ]  # fmt: skip
    by_type = (
        patients[patients.cancer_type.isin(TYPES)]
        .groupby("cancer_type")
        .tp53.agg(["mean", "size"])
        .sort_values("mean")
        .rename(index=tcga.CANCER_NAMES)
    )
    mo.md(
        f"Across all {len(patients):,} sequenced tumours: "
        f"**{100 * patients.tp53.mean():.1f}%** carry a TP53 mutation."
    )
    return (by_type,)


@app.cell
def _(COLOR, by_type, plt):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.barh(by_type.index, 100 * by_type["mean"], color=COLOR.BLUE.value,
            height=0.7)  # fmt: skip
    for _i, (_v, _n) in enumerate(zip(by_type["mean"], by_type["size"])):
        ax.text(100 * _v + 1, _i, f"{100 * _v:.0f}%  (n={_n})",
                va="center", fontsize=7)  # fmt: skip
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of tumours with a TP53 mutation")
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Which other genes are mutated in more than 10% of all tumours? Count
       patients per `gene_symbol`. Watch out for very long genes such as
       TTN: why do they top the list without being cancer genes?
    2. In breast cancer, compare TP53 mutation rates in tumours with and
       without a PIK3CA mutation.
    """)
    return


if __name__ == "__main__":
    app.run()
