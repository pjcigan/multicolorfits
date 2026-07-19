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
# # Getting started (scripting)
#
# Stretch, colorize, and combine FITS layers with the v3 session API.  For the
# full NGC 602 narrative (OpenFITS download, WCS rotation, inverse backgrounds)
# see [ngc602_same_grid](ngc602_same_grid.py).  Conceptual guides:
# [color compositing](../guide/color_compositing.md),
# [intensity scaling](../guide/intensity_scaling.md),
# [session model](../guide/session_model.md),
# [WCS & alignment](../guide/wcs_and_alignment.md).

# %%
import multicolorfits as mcf
from multicolorfits.figures import make_combined_figure

# %% [markdown]
# ## Minimal `McfSession` recipe
#
# Replace the paths with three registered FITS files on your machine.  Layers
# must share a pixel grid (use **Align layers…** in the GUI or `align_stack`
# when they do not).

# %%
# s = mcf.McfSession(n_panels=3)
# s.load_files(
#     ['layer_a.fits', 'layer_b.fits', 'layer_c.fits'],
#     colors=['#BE599E', '#DEA215', '#77C0F9'],
#     labels=['A', 'B', 'C'],
# )
# s.compose.combine_mode = 'lab'      # recommended alternative to classic RGB sum
# s.compose.combine_blend = 'screen'
# s.compose.gamma = 2.2
# s.compose.title = 'My target'
# s.compose.show_legend = True
# s.compose.show_combo_swatch = True  # overlap preview matches combine_mode
#
# rgb = s.render_combined()                 # full float64 for export
# # rgb = s.render_combined(preview=True)   # fast float32 path while iterating
# fig, ax = make_combined_figure(s, combined=rgb)
# fig.savefig('combined.png', dpi=150, bbox_inches='tight')

# %% [markdown]
# ## One-shot helper (arrays already aligned)

# %%
# rgb = mcf.combine_layers(
#     [data_a, data_b, data_c],
#     colors=['#BE599E', '#DEA215', '#77C0F9'],
#     colorspace='lab', blend='screen', gamma=2.2)

# %% [markdown]
# ## GUI and export
#
# ```python
# mcf.gui(session=s)                 # browser GUI (needs [web] or [all])
# # mcf.gui_qt(session=s)            # desktop GUI (needs [qt] or [all])
# s.save_state('session.json')
# print(s.export_script())           # recreate the plot without the GUI
# ```
#
# Need more than four layers?  Both GUIs ship with **+ Add panel** (up to 16),
# or call `s.add_panel()` in scripts.

# %%
print('multicolorfits', mcf.__version__)
