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
# # NGC 602 — same pixel grid
#
# Download IR / R / B OpenFITS products from
# [Chandra OpenFITS](http://chandra.harvard.edu/photo/openFITS/multiwavelength_data.html).
# The frames already share a grid.  Full narrative:
# [examples/ngc602](../examples/ngc602.md).

# %%
import multicolorfits as mcf
from multicolorfits.figures import make_combined_figure

FILES = [
    './ngc602_ir.fits',
    './ngc602_optical_R.fits',
    './ngc602_optical_B.fits',
]
COLORS = ['#BE599E', '#DEA215', '#77C0F9']  # purple / orange / blue
LABELS = ['IR', 'R', 'B']

# %% [markdown]
# ## Session + Lab compositing

# %%
# s = mcf.McfSession(n_panels=3)
# s.load_files(FILES, colors=COLORS, labels=LABELS)
# s.compose.combine_mode = 'lab'
# s.compose.combine_blend = 'screen'
# s.compose.gamma = 2.2
# s.compose.title = 'NGC 602'
# s.compose.show_legend = True
# s.compose.show_combo_swatch = True
#
# rgb = s.render_combined()
# fig, ax = make_combined_figure(s, combined=rgb)
# fig.savefig('n602_POB_lab.png', dpi=150, bbox_inches='tight')

# %% [markdown]
# White canvas without the legacy hex-invert + `inverse=True` path:
#
# ```python
# s.compose.combine_background = 'white'  # or 'transparent'
# ```
#
# Interactive finishing:
#
# ```python
# mcf.gui(session=s)
# print(s.export_script())
# ```

# %% [markdown]
# ## Classic low-level pipeline
#
# Still fully supported for v2-compatible scripts.

# %%
# import astropy.io.fits as pyfits
#
# ir, irhdr = pyfits.getdata(FILES[0], header=True)
# R, Rhdr = pyfits.getdata(FILES[1], header=True)
# B, Bhdr = pyfits.getdata(FILES[2], header=True)
#
# layers = []
# for data, color in zip([ir, R, B], COLORS):
#     grey = mcf.to_grey_rgb(data, rescalefn='linear')
#     layers.append(mcf.colorize_image(grey, color, colorintype='hex'))
#
# pob_rgb = mcf.combine_multicolor(layers, gamma=2.2)
# pob_lab = mcf.combine_multicolor_colorspace(
#     layers, colorspace='lab', blend='screen', gamma=2.2)
#
# mcf.plot_combined_rgb(
#     pob_lab, irhdr, 'NGC 602', 'n602_POB.jpg',
#     tickcolor='#D9D5C5', labelcolor='k', facecolor='w', minorticks=True)

# %% [markdown]
# ## North-up WCS
#
# OpenFITS NGC 602 is rotated ~90° from north.  Either edit a dummy header and
# `reproject_image`, or (with `[reproject]`) align to ICRS from the session / GUI:
#
# ```python
# s.align_panels(target='icrs')   # or GUI: Align layers… → ICRS
# ```
#
# Side-by-side Lab / RGB / HSV grids: [`examples/colorspace_comparison.md`](../../examples/colorspace_comparison.md).

# %%
print('Uncomment the cells above after downloading the OpenFITS products.')
