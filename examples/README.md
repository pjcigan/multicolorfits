# Examples

Worked narratives for Sphinx docs live under
[`docs/examples/`](../docs/examples/index.md) (NGC 602, WLM, M74
transparent cutouts, color-space suite, gallery). This directory keeps
**generators** and short GitHub pointers.

| File | Role |
|------|------|
| [`compare_colorspaces.py`](compare_colorspaces.py) | Builds light/dark PNG grids → `output/` |
| [`wide_gamut_probe.py`](wide_gamut_probe.py) | Research probe (not a user feature) |
| [`ngc602.md`](ngc602.md), [`wlm.md`](wlm.md), [`colorspace_comparison.md`](colorspace_comparison.md) | Pointers → Sphinx examples |

```bash
python examples/compare_colorspaces.py
python docs/make_docs_figures.py
python docs/make_showcase_figures.py   # real-data NGC 602 / M74 showcase figures
```

Set ``MCF_KEPLER_DIR`` / ``MCF_NGC602_DIR`` if your local FITS caches are not
in the default locations.
