# CLAUDE.md

Marimo notebooks reproducing the TCGA figures of the Cancer and the Cell Cycle lectures (Lectures 2 and 3). Public teaching repo: students clone it and run the notebooks.

## Figure Style — Non-Negotiable Defaults

**Always** apply the lab matplotlib style before creating any figure:

```python
import matplotlib.pyplot as plt

import tcga

tcga.style()  # lab style, scaled up for screen reading
```

**Always** use the `COLOR` enum from `colors.py` for plot colors:

```python
from colors import COLOR

ax.plot(x, y, color=COLOR.BLUE.value)
ax.scatter(x, y, color=COLOR.PINK.value)
```

Never use raw hex strings or default matplotlib colors in figure code. If a new color is needed, add it to the `COLOR` enum in `colors.py` first.

## Color Palette Reference

| Name | Hex | Typical use |
|------|-----|-------------|
| `COLOR.BLUE` | `#526C94` | Primary data series |
| `COLOR.LIGHT_BLUE` | `#75B1CE` | Secondary / paired condition |
| `COLOR.PINK` | `#DC6B83` | Contrast / highlight |
| `COLOR.YELLOW` | `#D8C367` | Third condition |
| `COLOR.TURQUOISE` | `#00bfb2` | Fourth condition |
| `COLOR.LIGHT_GREEN` | `#CCDBA2` | Fifth condition |
| `COLOR.LAVENDER` | `#C6B2D1` | Sixth condition |
| `COLOR.PURPLE` | `#654875` | Seventh condition |
| `COLOR.OLIVE` | `#889466` | Eighth condition |
| `COLOR.GREY` | `#D4D3CF` | Background / negative control |
| `COLOR.DARKGREY` | `#4A4A4A` | Text / annotations |

## Data Analysis Conventions

- **DataFrames**: Use `pandas` for tabular data; name frames descriptively (`cell_df`, `intensity_df`)
- **Plots**: `matplotlib` for fine-grained control, `seaborn` for statistical plots — always pass `palette` from `COLOR` values
- **Stats**: State the test used and its assumptions before reporting p-values
- **Images**: `scikit-image` or `tifffile` for microscopy images (add if needed)
- **Units**: Always label axes with units; use µm not um, s not sec

## Project Layout

```
ccc-tcga-figures/
├── hhlab_style01.mplstyle   # matplotlib style — never edit
├── notebooks/
│   ├── NN_<topic>.py         # one marimo notebook per figure
│   ├── tcga.py               # loading, slide definitions, style()
│   └── colors.py             # COLOR enum — extend here for new colors
├── scripts/build_extract.py  # builds data/ from tcga-db (~/code/tcga-db)
├── data/                     # committed parquet extract (~16 MB) — keep small
└── results/                  # figures and tables written by notebooks
```

## Conventions

- Numbers must match the lecture slides: definitions mirror
  `~/webdev/cancer-cell-cycle/scripts/tcga_export*.py`.
- marimo: cell-local names start with `_`; last expression is the output.
- Check: `uv run ruff check && uv run marimo check notebooks/[0-9]*.py`.

## Environment

- `uv add <package>` to add dependencies — never `pip install`
- `uv run marimo edit notebooks/<name>.py`
