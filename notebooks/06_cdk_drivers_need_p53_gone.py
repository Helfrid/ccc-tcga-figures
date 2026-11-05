import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 6. Tumours driving CDK2 have lost p53

    **Lecture 3.** Deregulated CDK2 (cyclin E amplified, MYC driving E2F)
    causes replication stress; p53 stops such cells. The model predicts
    that tumours with these drivers must have lost p53.

    For each amplified gene: the fraction of tumours with TP53 mutated
    (**observed**), against the fraction **expected** from their cancer
    types alone (each tumour weighted by the TP53 rate in its type). This
    controls for cancer types that simply have many TP53 mutations.

    MDM2 is the control the other way round: it destroys p53, so tumours
    that amplify MDM2 rarely need a TP53 mutation as well.
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
    _m = tcga.load("mutations")
    _coding = tcga.MISSENSE + tcga.TRUNCATING
    _tp53 = set(
        _m.patient_id[
            (_m.gene_symbol == "TP53")
            & _m.variant_classification.isin(_coding)
        ]
    )
    cn = tcga.copy_number_wide()
    tumours = tcga.sequenced().set_index("patient_id")
    tumours = tumours[tumours.index.isin(cn.index)]
    tumours["tp53"] = tumours.index.isin(_tp53)
    # the TP53 rate of each tumour's own cancer type
    tumours["type_rate"] = tumours.groupby("cancer_type").tp53.transform(
        "mean"
    )
    tumours = tumours.join(cn)
    return (tumours,)


@app.cell
def _(mo):
    threshold = mo.ui.slider(
        0.5, 2.0, step=0.1, value=1.0, show_value=True,
        label="amplified: log2 copy ratio above",
    )  # fmt: skip
    threshold
    return (threshold,)


@app.cell
def _(pd, threshold, tumours):
    GENES = {"CCNE1": "cyclin E", "MYC": "MYC", "CCND1": "cyclin D1",
             "MDM2": "MDM2"}  # fmt: skip
    _rows = []
    for _g, _label in GENES.items():
        _amp = tumours[tumours[_g] > threshold.value]
        _rows.append({
            "gene": _label if _label == _g else f"{_label} ({_g})",
            "n": len(_amp),
            "observed": _amp.tp53.mean(), "expected": _amp.type_rate.mean(),
        })  # fmt: skip
    table = pd.DataFrame(_rows).set_index("gene")
    table.round(3)
    return (table,)


@app.cell
def _(COLOR, np, plt, table):
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    _y = np.arange(len(table))[::-1]
    ax.barh(_y, 100 * table.observed, height=0.6, color=COLOR.BLUE.value,
            label="observed")  # fmt: skip
    ax.scatter(100 * table.expected, _y, marker="|", s=500, color="black",
               linewidths=2, label="expected from cancer type")  # fmt: skip
    ax.set_yticks(_y, [f"{g}  n={n}" for g, n in table.n.items()])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% with TP53 mutated, among tumours with the gene amplified")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Move the slider. Does the pattern hold for stronger amplifications?
    2. Add CDK4 and MDM4. Which behaves like MDM2, and why?
    3. Co-occurrence is not order. What experiment would show that cyclin E
       comes first and p53 loss second?
    """)
    return


if __name__ == "__main__":
    app.run()
