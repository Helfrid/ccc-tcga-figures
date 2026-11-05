# /// script
# requires-python = ">=3.12"
# dependencies = ["duckdb>=1.1"]
# ///
"""Build the small parquet extract in data/ from a full tcga-db database.

You do not need to run this: the extract is committed. It is here so every
number can be traced back to the full database
(https://github.com/Helfrid/tcga-db).

    uv run scripts/build_extract.py [path/to/tcga.duckdb]

Everything is per patient. A few patients have more than one tumour sample;
for copy number, ploidy and arm calls we keep the highest value, as on the
slides.

Files written:
    patients.parquet        patient_id, cancer type (COAD + READ merged into
                            COADREAD), sequenced (has mutation calls)
    mutations.parquet       protein-changing mutations (MC3 public MAF)
    copy_number.parquet     gene-level log2 copy ratio for ~60 cancer and
                            cell cycle genes (Broad SNP6, gene level)
    genome_altered.parquet  fraction of autosomal genes (every 5th, ~3,600)
                            with |log2 copy ratio| > 0.3
    arms.parquet            median log2 copy ratio of each chromosome arm
                            (every 3rd gene; hg38 centromeres)
    ploidy.parquet          ABSOLUTE ploidy, genome doublings, purity
    genes.parquet           protein-coding genes with chromosome arm
    expression_brca.parquet breast tumours, log2(TPM + 0.001) for all 8q genes
                            and 1,500 random genes elsewhere
"""

from __future__ import annotations

import logging
import random
import sys
from pathlib import Path

import duckdb

logger = logging.getLogger("build_extract")

DB_DEFAULT = Path.home() / "code/tcga-db/data/tcga.duckdb"
OUT = Path(__file__).resolve().parent.parent / "data"

CODING = (
    "Missense_Mutation",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
)

CN_GENES = [
    # Lecture 2: growth signals and the Rb pathway
    "EGFR",
    "ERBB2",
    "KRAS",
    "NRAS",
    "HRAS",
    "BRAF",
    "MYC",
    "MYCN",
    "MYCL",
    "CCND1",
    "CCND2",
    "CCND3",
    "CDK4",
    "CDK6",
    "CDKN2A",
    "CDKN2B",
    "RB1",
    "E2F1",
    "E2F3",
    "PIK3CA",
    "PTEN",
    "AKT1",
    "MET",
    "FGFR1",
    "FGFR2",
    "KIT",
    "PDGFRA",
    "NF1",
    "APC",
    "SMAD4",
    # Lecture 3: p53, checkpoints, mitosis
    "TP53",
    "MDM2",
    "MDM4",
    "CDKN1A",
    "CDKN1B",
    "CCNE1",
    "CCNE2",
    "CDK2",
    "CCNA2",
    "CCNB1",
    "CDK1",
    "ATR",
    "ATM",
    "CHEK1",
    "CHEK2",
    "WEE1",
    "PKMYT1",
    "BRCA1",
    "BRCA2",
    "PARP1",
    "AURKA",
    "AURKB",
    "PLK1",
    "MASTL",
    "KIF18A",
    "MAD2L1",
    "BUB1B",
    "TERT",
]

# hg38 centromere midpoints (Mb); p arms of acrocentric chromosomes left out
CENTROMERE: dict[int, float] = {
    1: 123.4, 2: 93.9, 3: 90.9, 4: 50.0, 5: 48.8, 6: 59.8, 7: 60.1, 8: 45.2,
    9: 43.0, 10: 39.8, 11: 53.4, 12: 35.5, 13: 17.7, 14: 17.2, 15: 19.0,
    16: 36.8, 17: 25.1, 18: 18.5, 19: 26.2, 20: 28.1, 21: 12.0, 22: 15.0,
}  # fmt: skip
ACROCENTRIC = {13, 14, 15, 21, 22}
ALTERED_LOG2 = 0.3


def arm_of(chrom: str | None, start: int | None) -> str | None:
    """Chromosome arm ('8q') of a gene, or None off the autosomes."""
    if not chrom or start is None:
        return None
    c = chrom.removeprefix("chr")
    if not c.isdigit():
        return None
    n = int(c)
    arm = "p" if start / 1e6 < CENTROMERE[n] else "q"
    return None if n in ACROCENTRIC and arm == "p" else f"{n}{arm}"


def copy(con: duckdb.DuckDBPyConnection, sql: str, name: str) -> None:
    """Write a query to data/<name>.parquet."""
    path = OUT / f"{name}.parquet"
    con.execute(f"copy ({sql}) to '{path}' (format parquet, compression zstd)")
    n = con.execute(f"select count(*) from '{path}'").fetchone()[0]
    logger.info(
        "%-24s %9d rows %7.1f MB", path.name, n, path.stat().st_size / 1e6
    )


def main() -> None:
    """Query tcga-db and write the parquet extract."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    db = Path(sys.argv[1]) if len(sys.argv) > 1 else DB_DEFAULT
    con = duckdb.connect(str(db), read_only=True)
    con.execute("set preserve_insertion_order = false")
    OUT.mkdir(exist_ok=True)

    con.execute(
        """
        create temp view tumour as
        select s.sample_id, s.patient_id from samples s where s.is_tumour
        """
    )
    copy(
        con,
        """
        select p.patient_id,
               case when cancer_type in ('COAD', 'READ') then 'COADREAD'
                    else cancer_type end as cancer_type,
               p.patient_id in (select patient_id from somatic_mutations)
                   as sequenced
        from patients p where cancer_type <> 'CONTROL' order by 1
        """,
        "patients",
    )
    copy(
        con,
        f"""
        select distinct patient_id, gene_symbol, variant_classification,
               hgvsp_short as protein_change
        from somatic_mutations where variant_classification in {CODING}
        order by 2, 1
        """,
        "mutations",
    )

    genes = con.execute(
        "select gene_id, gene_symbol, chromosome, chrom_start, "
        "is_protein_coding from genes"
    ).fetchall()
    cn_cols = [
        r[0] for r in con.execute("describe copy_number_wide").fetchall()
    ][1:]
    sym = {g: s for g, s, *_ in genes}
    pos = {g: (c, st) for g, _, c, st, _ in genes}

    wanted = {s: g for g, s in sym.items() if s in CN_GENES and g in cn_cols}
    missing = set(CN_GENES) - set(wanted)
    if missing:
        logger.warning("no copy number for %s", sorted(missing))
    sel = ", ".join(f"{g} as {s}" for s, g in wanted.items())
    names = ", ".join(wanted)
    copy(
        con,
        f"""
        select patient_id, gene, round(max(log2), 3) as log2_ratio
        from (unpivot (select t.patient_id, {sel}
                       from copy_number_wide join tumour t using (sample_id))
              on {names} into name gene value log2)
        group by 1, 2 order by 2, 1
        """,
        "copy_number",
    )

    auto = [g for g in cn_cols if arm_of(*pos.get(g, (None, None)))][::5]
    cols = ", ".join(auto)
    copy(
        con,
        f"""
        select patient_id, round(max(f), 4) as fraction_altered
        from (
            select sample_id, avg((abs(v) > {ALTERED_LOG2})::int) as f
            from (unpivot (select sample_id, {cols} from copy_number_wide)
                  on {cols} into name gene value v)
            group by 1
        ) join tumour using (sample_id) group by 1 order by 1
        """,
        "genome_altered",
    )

    by_arm: dict[str, list[str]] = {}
    for g in cn_cols:
        a = arm_of(*pos.get(g, (None, None)))
        if a:
            by_arm.setdefault(a, []).append(g)
    # one row-wise median per arm, then long format
    meds = ", ".join(
        f"list_median([{', '.join(gs[::3])}]) as a{a}"
        for a, gs in by_arm.items()
    )
    arms = ", ".join(f"a{a}" for a in by_arm)
    copy(
        con,
        f"""
        select patient_id, substr(arm, 2) as arm,
               round(max(v), 3) as log2_ratio
        from (unpivot (select sample_id, {meds} from copy_number_wide)
              on {arms} into name arm value v)
        join tumour using (sample_id) group by 1, 2 order by 1, 2
        """,
        "arms",
    )

    copy(
        con,
        """
        select t.patient_id, max(q.ploidy) as ploidy,
               max(q.genome_doublings) as genome_doublings,
               max(q.tumour_purity) as purity
        from sample_qc q join tumour t using (sample_id)
        where q.genome_doublings is not null group by 1 order by 1
        """,
        "ploidy",
    )

    rows = [
        (g, s, c, st, arm_of(c, st))
        for g, s, c, st, pc in genes
        if pc and arm_of(c, st)
    ]
    con.execute(
        "create temp table pc (gene_id varchar, gene_symbol varchar, "
        "chromosome varchar, start bigint, arm varchar)"
    )
    con.executemany("insert into pc values (?, ?, ?, ?, ?)", rows)
    copy(con, "select * from pc order by chromosome, start", "genes")

    ex = {r[0] for r in con.execute("describe expression_tpm_wide").fetchall()}
    on8q = [g for g, _, _, _, a in rows if a == "8q" and g in ex]
    rest = [g for g, _, _, _, a in rows if a not in ("8p", "8q") and g in ex]
    random.seed(2)
    keep = on8q + random.sample(rest, 1500)
    sel = ", ".join(f"round({g}, 3) as {g}" for g in keep)
    copy(
        con,
        f"""
        select t.patient_id, {sel}
        from expression_tpm_wide join tumour t using (sample_id)
        join patients p using (patient_id) where p.cancer_type = 'BRCA'
        """,
        "expression_brca",
    )


if __name__ == "__main__":
    main()
