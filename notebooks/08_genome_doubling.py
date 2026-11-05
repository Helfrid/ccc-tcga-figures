import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 8. Doubled genomes

    **Lecture 3.** ABSOLUTE estimates each tumour's ploidy (average number
    of copies of the genome) and how many times the whole genome has been
    doubled. A cell that fails cytokinesis, or skips mitosis, becomes
    tetraploid (4n). Healthy cells stop after that, through p53; tumours
    that have lost p53 carry on.

    Three questions:

    1. How common is doubling in each cancer type?
    2. Is it more common in TP53-mutant tumours, *within* each type?
    3. Doubled tumours started at 4n. Where are they now?
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
    ploidy = tcga.load("ploidy").merge(tcga.load("patients"))
    ploidy["doubled"] = ploidy.genome_doublings > 0
    ploidy.head()
    return (ploidy,)


@app.cell
def _(COLOR, ploidy, plt, tcga):
    TYPES = [
        "OV", "ESCA", "LUSC", "HNSC", "PAAD", "COADREAD", "STAD", "BLCA",
        "LUAD", "BRCA", "UCEC", "LIHC", "GBM", "SARC", "SKCM", "PRAD", "KIRC",
        "THCA",
    ]  # fmt: skip
    _w = (
        ploidy[ploidy.cancer_type.isin(TYPES)]
        .groupby("cancer_type")
        .doubled.mean()
        .sort_values()
        .rename(index=tcga.CANCER_NAMES)
    )
    fig1, _ax = plt.subplots(figsize=(6, 5))
    _ax.barh(_w.index, 100 * _w, color=COLOR.PURPLE.value, height=0.7)
    _ax.set_xlim(0, 100)
    _ax.set_xlabel("% of tumours with a doubled genome")
    fig1.tight_layout()
    fig1
    return


@app.cell
def _(mo, ploidy, tcga):
    _m = tcga.load("mutations")
    _coding = tcga.MISSENSE + tcga.TRUNCATING
    _tp53 = set(
        _m.patient_id[
            (_m.gene_symbol == "TP53")
            & _m.variant_classification.isin(_coding)
        ]
    )
    seq = ploidy[ploidy.patient_id.isin(tcga.sequenced().patient_id)].copy()
    seq["tp53"] = seq.patient_id.isin(_tp53)

    _pan = seq.groupby("tp53").doubled.mean()
    _g = seq.groupby(["cancer_type", "tp53"]).doubled
    _f, _n = _g.mean().unstack(), _g.size().unstack()
    _ok = (_n[True] >= 20) & (_n[False] >= 20)
    _higher = (_f[True] > _f[False])[_ok].sum()
    mo.md(
        f"Doubled: **{100 * _pan[True]:.0f}%** of TP53-mutant tumours against "
        f"**{100 * _pan[False]:.0f}%** of wild-type ones. Within cancer types "
        f"(at least 20 per group) higher in the mutant in **{_higher} of "
        f"{_ok.sum()}**, so tissue mix does not explain it."
    )
    return


@app.cell
def _(COLOR, mo, np, ploidy, plt):
    _once = ploidy[ploidy.genome_doublings == 1].ploidy
    _bins = np.arange(1, 6.01, 0.1)
    fig2, ax2 = plt.subplots(figsize=(8, 3.6))
    ax2.hist(ploidy.ploidy[~ploidy.doubled], bins=_bins,
             color=COLOR.LIGHT_BLUE.value, label="not doubled")  # fmt: skip
    ax2.hist(ploidy.ploidy[ploidy.doubled], bins=_bins,
             color=COLOR.PURPLE.value, alpha=0.8, label="doubled")  # fmt: skip
    ax2.axvline(4, color="black", ls="--", lw=1)
    ax2.text(4.05, ax2.get_ylim()[1] * 0.9, "4n at the moment of doubling")
    ax2.axvline(_once.median(), color=COLOR.PINK.value, lw=1.5)
    ax2.text(_once.median() - 0.05, ax2.get_ylim()[1] * 0.75,
             f"once doubled, now: {_once.median():.2f} (median)",
             ha="right", color=COLOR.PINK.value)  # fmt: skip
    ax2.set_xlabel("tumour ploidy (average copies of the genome)")
    ax2.set_ylabel("tumours")
    ax2.legend(loc="upper right")
    fig2.tight_layout()
    mo.vstack([
        fig2,
        mo.md(
            f"{len(ploidy):,} tumours, {100 * ploidy.doubled.mean():.0f}% "
            f"doubled. Of the once-doubled tumours only "
            f"**{100 * (_once >= 4).mean():.0f}%** are still at 4 or above: "
            "after doubling, chromosomes are lost."
        ),
    ])  # fmt: skip
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Is a doubled genome bad news for the patient? Survival data are not
       in this extract, but they are in tcga-db (`survival` table): try it
       there, within one cancer type.
    2. Do twice-doubled tumours (`genome_doublings == 2`) sit near 8n?
    """)
    return


if __name__ == "__main__":
    app.run()
