import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 1. Which oncogenes drive which cancers?

    **Lecture 2.** For each cancer type: the percent of tumours with each
    oncogene switched on, by an activating mutation or by amplification.

    Definitions (the same as on the slides):

    - **mutated**: a missense or in-frame mutation at a known activating
      codon (KRAS G12, BRAF V600, the EGFR kinase domain, ...);
    - **amplified**: gene-level log2 copy ratio above 1 (about four copies
      or more). MYC, CCND1 and CDK4 count only when amplified;
    - denominator: tumours with both mutation calls and copy number.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    import tcga

    tcga.style()
    return mo, np, pd, plt, tcga


@app.cell
def _(tcga):
    mutations = tcga.load("mutations")
    cn = tcga.copy_number_wide()
    patients = tcga.sequenced()
    # tumours with both kinds of data
    patients = patients[patients.patient_id.isin(cn.index)]
    patients.cancer_type.value_counts().head()
    return cn, mutations, patients


@app.cell
def _(cn, mutations, tcga):
    GENES = ["EGFR", "ERBB2", "KRAS", "NRAS", "BRAF", "MYC", "CCND1", "CDK4"]
    AMPLIFIED_ONLY = {"MYC", "CCND1", "CDK4"}
    BY_AMPLIFICATION = {"EGFR", "ERBB2", "MYC", "CCND1", "CDK4"}

    altered = {}
    for _g in GENES:
        _mut = set() if _g in AMPLIFIED_ONLY else tcga.mutated(mutations, _g)
        _amp = tcga.amplified(cn, _g) if _g in BY_AMPLIFICATION else set()
        altered[_g] = _mut | _amp
    {g: len(p) for g, p in altered.items()}
    return GENES, altered


@app.cell
def _(GENES, altered, patients, pd, tcga):
    TYPES = [
        "LUAD", "LUSC", "PAAD", "COADREAD", "STAD", "ESCA", "SKCM",
        "THCA", "BRCA", "HNSC", "BLCA", "UCEC", "GBM", "SARC",
    ]  # fmt: skip
    _rows = {}
    for _t in TYPES:
        _pts = set(patients.patient_id[patients.cancer_type == _t])
        _rows[tcga.CANCER_NAMES[_t]] = {
            _g: 100 * len(altered[_g] & _pts) / len(_pts) for _g in GENES
        } | {"n": len(_pts)}
    landscape = pd.DataFrame(_rows).T
    landscape.round(1)
    return (landscape,)


@app.cell
def _(GENES, landscape, np, plt):
    _v = landscape[GENES]
    fig, ax = plt.subplots(figsize=(7, 5.5))
    _im = ax.imshow(_v, cmap="Blues", vmin=0, vmax=70, aspect="auto")
    ax.set_xticks(range(len(GENES)), GENES)
    ax.set_yticks(
        range(len(_v)), [f"{t} ({n:.0f})" for t, n in landscape.n.items()]
    )
    ax.xaxis.tick_top()
    for (_i, _j), _x in np.ndenumerate(_v.to_numpy()):
        if _x >= 1:
            ax.text(
                _j, _i, f"{_x:.0f}", ha="center", va="center", fontsize=7,
                color="white" if _x > 40 else "black",
            )  # fmt: skip
    fig.colorbar(_im, ax=ax, shrink=0.6, label="% of tumours altered")
    ax.spines[:].set_visible(False)
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. KRAS is mutated in most pancreatic cancers. Which codon is hit most
       often? (Filter `mutations` for KRAS and count `protein_change`.)
    2. Change the amplification threshold `tcga.AMP_LOG2` to 0.5. Which
       genes change most, and why might a lower cut-off be misleading?
    3. Add PIK3CA to the heatmap. Which codons would you call activating?
    """)
    return


if __name__ == "__main__":
    app.run()
