"""Shared data loading, definitions and plot style for the notebooks.

Every notebook starts with ``import tcga``. The definitions here are the
ones used on the lecture slides, so the numbers match.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# --- definitions used on the slides ---------------------------------------

AMP_LOG2 = 1.0
"""Amplified: gene-level log2 copy ratio above 1 (about 4+ copies)."""

DEL_LOG2 = -1.0
"""Deleted: log2 copy ratio below -1 (both copies lost in most cells)."""

MISSENSE = ("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")
TRUNCATING = (
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
)

HOTSPOTS: dict[str, set[int]] = {
    "KRAS": {12, 13, 61, 117, 146},
    "NRAS": {12, 13, 61, 117, 146},
    "HRAS": {12, 13, 61},
    "BRAF": {464, 466, 469, 594, 597, 600, 601},
    # kinase domain (exons 18-21) plus the glioma extracellular hotspots
    "EGFR": set(range(688, 876)) | {108, 289, 596, 598},
    "ERBB2": {310, 755, 777, 842},
}
"""Activating codons counted as 'mutated' for each oncogene."""

CANCER_NAMES: dict[str, str] = {
    "ACC": "Adrenocortical",
    "BLCA": "Bladder",
    "BRCA": "Breast",
    "CESC": "Cervix",
    "CHOL": "Bile duct",
    "COADREAD": "Colorectal",
    "DLBC": "Lymphoma (DLBC)",
    "ESCA": "Oesophagus",
    "GBM": "Glioblastoma",
    "HNSC": "Head and neck",
    "KICH": "Kidney chromophobe",
    "KIRC": "Kidney",
    "KIRP": "Kidney papillary",
    "LAML": "Leukaemia (AML)",
    "LGG": "Lower-grade glioma",
    "LIHC": "Liver",
    "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous",
    "MESO": "Mesothelioma",
    "OV": "Ovary",
    "PAAD": "Pancreas",
    "PCPG": "Phaeochromocytoma",
    "PRAD": "Prostate",
    "SARC": "Sarcoma",
    "SKCM": "Melanoma",
    "STAD": "Stomach",
    "TGCT": "Testis",
    "THCA": "Thyroid",
    "THYM": "Thymus",
    "UCEC": "Endometrium",
    "UCS": "Uterine carcinosarcoma",
    "UVM": "Uveal melanoma",
}

SOURCE = (
    "TCGA PanCancer Atlas via UCSC Xena (tcga-db): MC3 mutations, "
    "Broad SNP6 copy number, ABSOLUTE ploidy, Toil RNA-seq."
)

# --- loading ----------------------------------------------------------------


def load(name: str) -> pd.DataFrame:
    """Read one table of the extract, e.g. ``load("mutations")``."""
    return pd.read_parquet(DATA / f"{name}.parquet")


def sequenced() -> pd.DataFrame:
    """Patients with mutation calls. Ovarian cancer is left out: its
    mutation calls cover too few patients to use."""
    p = load("patients")
    return p[p.sequenced & (p.cancer_type != "OV")]


def copy_number_wide() -> pd.DataFrame:
    """Copy number as one row per patient, one column per gene."""
    return load("copy_number").pivot(
        index="patient_id", columns="gene", values="log2_ratio"
    )


# --- calling alterations ------------------------------------------------------


def codon(change: str | None) -> int | None:
    """First residue number in a protein change: 'p.G12D' -> 12."""
    if not change:
        return None
    m = re.search(r"(\d+)", change)
    return int(m.group(1)) if m else None


def mutated(mutations: pd.DataFrame, gene: str) -> set[str]:
    """Patients with an activating (hotspot) mutation in an oncogene."""
    m = mutations[
        (mutations.gene_symbol == gene)
        & mutations.variant_classification.isin(MISSENSE)
    ]
    hits = m.protein_change.map(codon).isin(HOTSPOTS.get(gene, set()))
    return set(m.patient_id[hits])


def inactivated(mutations: pd.DataFrame, gene: str) -> set[str]:
    """Patients with a mutation that inactivates a tumour suppressor:
    truncating, or any coding change for CDKN2A (missense breaks p16)."""
    kinds = TRUNCATING + (MISSENSE if gene == "CDKN2A" else ())
    m = mutations[
        (mutations.gene_symbol == gene)
        & mutations.variant_classification.isin(kinds)
    ]
    return set(m.patient_id)


def amplified(cn: pd.DataFrame, gene: str) -> set[str]:
    """Patients with the gene amplified (``cn`` from copy_number_wide)."""
    return set(cn.index[cn[gene] > AMP_LOG2])


def deleted(cn: pd.DataFrame, gene: str) -> set[str]:
    """Patients with the gene deleted on both copies."""
    return set(cn.index[cn[gene] < DEL_LOG2])


# --- style --------------------------------------------------------------------


def style() -> None:
    """Lab figure style, scaled up for reading on screen."""
    plt.style.use(ROOT / "hhlab_style01.mplstyle")
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.title_fontsize": 8,
            "legend.frameon": False,
            "legend.scatterpoints": 1,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.5,
            "figure.subplot.left": 0.125,
            "figure.subplot.right": 0.95,
            "figure.subplot.bottom": 0.12,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        }
    )
