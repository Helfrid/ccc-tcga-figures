import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # 3. Where are genes hit? Oncogenes and tumour suppressors

    **Lectures 2 and 3.** A *lollipop plot*: the protein runs left to right;
    each stem is one residue, its height the number of patients with a
    mutation there. Missense changes point up, truncating ones (stop,
    frameshift, splice) point down.

    - An **oncogene** (KRAS, BRAF) is hit at a few exact spots, always
      missense: only those changes switch the protein on.
    - A **tumour suppressor** (RB1) is broken anywhere, mostly truncated:
      any change that destroys it will do.
    - **TP53** is the odd one: a suppressor hit at hotspots (R175, R248,
      R273), because the mutant protein also poisons the normal copy.
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
def _(mo):
    PROTEINS = {
        "KRAS": {"length": 189, "domains": [(1, 166, "GTPase")]},
        "BRAF": {"length": 766,
                 "domains": [(155, 227, "RAS-binding"), (457, 717, "Kinase")]},
        "RB1": {"length": 928, "domains": [(373, 771, "Pocket")]},
        "TP53": {"length": 393,
                 "domains": [(1, 61, "Transactivation"),
                             (102, 292, "DNA-binding"),
                             (325, 356, "Tetramer")]},
    }  # fmt: skip
    gene = mo.ui.dropdown(list(PROTEINS), value="TP53", label="Gene")
    gene
    return PROTEINS, gene


@app.cell
def _(PROTEINS, gene, tcga):
    _m = tcga.load("mutations")
    _m = _m[_m.gene_symbol == gene.value].copy()
    _m["pos"] = _m.protein_change.map(tcga.codon)
    _m["kind"] = _m.variant_classification.map(
        lambda v: "missense" if v in tcga.MISSENSE
        else "truncating" if v in tcga.TRUNCATING else None
    )  # fmt: skip
    _m = _m.dropna(subset=["pos", "kind"])
    _m = _m[_m.pos <= PROTEINS[gene.value]["length"]]
    # patients per residue and kind
    sites = (
        _m.groupby(["pos", "kind"]).patient_id.nunique().unstack(fill_value=0)
    )
    sites.sort_values("missense", ascending=False).head(8)
    return (sites,)


@app.cell
def _(COLOR, PROTEINS, gene, plt, sites):
    _spec = PROTEINS[gene.value]
    fig, ax = plt.subplots(figsize=(9, 3.4))
    for _kind, _sign, _c in [("missense", 1, COLOR.BLUE.value),
                             ("truncating", -1, COLOR.PINK.value)]:  # fmt: skip
        if _kind not in sites:
            continue
        _s = sites[_kind][sites[_kind] > 0]
        ax.vlines(_s.index, 0, _sign * _s, color=_c, lw=1)
        ax.scatter(_s.index, _sign * _s, s=10, color=_c, label=_kind, zorder=3)
    # protein track below the truncating stems
    _h = sites.to_numpy().max() * 0.07
    _low = sites.get("truncating", sites.iloc[:, 0] * 0).max()
    _y = -_low - 2.5 * _h
    ax.add_patch(plt.Rectangle((1, _y - _h / 3), _spec["length"], _h / 1.5,
                               color=COLOR.GREY.value))  # fmt: skip
    for _a, _b, _name in _spec["domains"]:
        ax.add_patch(plt.Rectangle((_a, _y - _h / 2), _b - _a, _h,
                                   color=COLOR.YELLOW.value))  # fmt: skip
        ax.text((_a + _b) / 2, _y - _h, _name, ha="center", va="top",
                fontsize=8)  # fmt: skip
    # label the three biggest hotspots
    if "missense" in sites:
        for _pos, _n in sites.missense.nlargest(3).items():
            ax.annotate(f"{_pos:.0f}", (_pos, _n), xytext=(0, 4),
                        textcoords="offset points", ha="center")  # fmt: skip
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlim(0, _spec["length"] + 5)
    ax.set_xlabel(f"{gene.value} residue")
    ax.set_ylabel("patients (truncating below)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Try this

    1. Add PIK3CA (1,068 residues) or PTEN (403) to `PROTEINS`. Which looks
       like an oncogene and which like a suppressor?
    2. For TP53, what fraction of all mutations fall in the DNA-binding
       domain (residues 102 to 292)?
    """)
    return


if __name__ == "__main__":
    app.run()
