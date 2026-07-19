# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     kernelspec:
#       display_name: Python 3
#       name: python3
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
# ---

# %% [markdown]
# # Component mosaic
#
# Hero combined image + per-layer colorized strip on a shared WCS.  Guide
# (with light/dark layout figures):
# [figures_and_mosaics](../guide/figures_and_mosaics.md).

# %%
import multicolorfits as mcf

# %% [markdown]
# ## Basic call
#
# ```python
# s = mcf.McfSession()
# s.load_files([...], colors=[...], labels=[...])
# s.compose.combine_mode = 'lab'
#
# fig, axes = mcf.make_component_mosaic(
#     s,
#     components='top',      # strip above hero; also bottom / left / right
#     max_per_line=3,
#     ticks='plain',          # plain | minimal | full
#     overlays='hero',        # beam/compass/scale bar on hero only
# )
# fig.savefig('mosaic.png', dpi=150, bbox_inches='tight')
# ```
#
# `axes['combined']` is the hero; `axes['components']` is the list of layer
# panels.  Empty leftover cells in the strip stay blank for custom annotations.
#
# Local regen scratch: `examples/output/component_mosaic/` → copy into
# `docs/_static/mosaics/` via `docs/make_docs_figures.py` (not committed).

# %%
print('mcf.make_component_mosaic — see guide/figures_and_mosaics.md')
