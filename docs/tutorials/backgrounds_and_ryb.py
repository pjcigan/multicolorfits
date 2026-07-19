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
# # Backgrounds and RYB paint mixing
#
# White / transparent canvases and subtractive **RYB** mixing.  Concepts and
# figures: [color compositing](../guide/color_compositing.md).  CMYK is
# API-only in the combine-mode dropdown pending GUI polishing — available as
# `combine_mode='cmyk'`.
#
# ![RYB paint mix vs additive RGB](../_static/compositing/paintmix_demo_light.png)

# %%
import multicolorfits as mcf

# %% [markdown]
# ## Backgrounds (preferred over `inverse=True`)
#
# Classic `combine_multicolor(..., inverse=True)` plus `hex_complement()` was the v2
# route to a light canvas.  Prefer an explicit background:
#
# ```python
# s.compose.combine_background = 'white'
# s.compose.combine_background = 'transparent'
# s.compose.combine_background = '#F5F0E6'   # custom
# s.compose.facecolor = 'none'              # transparent figure canvas when saving
# ```
#
# In the browser GUI the Compositing tab exposes the same choices (and a
# side-by-side montage of white / RYB / transparent):
# [gui walkthrough](gui_web_walkthrough.md).

# %% [markdown]
# ## RYB (paint-like subtractive)
#
# ```python
# s.compose.combine_mode = 'ryb'
# # blend is ignored for RYB/CMYK — chroma mixes subtractively
# rgb = s.render_combined()
# ```
#
# Or with the low-level helpers:
#
# ```python
# rgb = mcf.combine_colorized_layers(colorized_list, mode='ryb', gamma=2.2)
# ```
#
# Yellow + blue → green on the paint wheel (not cyan-ish white as in additive
# RGB).  That is intentional.  For additive Lab mixing that keeps hues distinct
# under bright overlaps, use `combine_mode='lab'` + `combine_blend='screen'`.

# %% [markdown]
# ## When to use what
#
# | Goal | Setting |
# |------|---------|
# | General science color composite | Lab / screen |
# | Classic channel sum | RGB (`combine_multicolor`) |
# | Paint / ink aesthetic | RYB |
# | Print-ink invert feel | CMYK (API) |
# | Slides with irregular outline | transparent background + cutout export |

# %%
print('Backgrounds + RYB — see docs/guide/color_compositing.md')
