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
# # Save and export
#
# How to leave a session with rasters, RGB FITS (with provenance), transparent
# cutouts, session JSON, and a recreation script.  Related:
# [saving_outputs guide](../guide/saving_outputs.md),
# [component mosaic](component_mosaic.py).

# %%
import multicolorfits as mcf
from multicolorfits.figures import make_combined_figure

# %% [markdown]
# ## Assume a loaded session
#
# ```python
# s = mcf.McfSession()
# s.load_files(['a.fits', 'b.fits'], colors=['#C11B17', '#4CC417'], labels=['Hα', '[OIII]'])
# s.compose.combine_mode = 'lab'
# s.compose.title = 'My nebula'
# ```

# %% [markdown]
# ## Raster (format from extension)

# %%
# rgb = s.render_combined()
# fig, ax = make_combined_figure(s, combined=rgb)
# fig.savefig('out.png', dpi=300, bbox_inches='tight')   # also .jpg / .pdf / .svg …
#
# # Publication / stamp-friendly (no ticks):
# s.apply_bare_plot_style(transparent=True)
# fig2, ax2 = make_combined_figure(s)
# fig2.savefig('stamp.png', dpi=200, transparent=True, bbox_inches='tight')

# %% [markdown]
# ## RGB FITS + provenance

# %%
# s.save_rgb_fits('combined_rgb.fits')   # 3-plane cube + HISTORY cards
# # Low-level: mcf.save_rgb_fits(path, rgb_array, header)

# %% [markdown]
# ## Transparent cutout (slides / stamps)

# %%
# rgba = s.export_transparent_cutout(size=400, alpha_lo=55, alpha_hi=99.5)
# # Or GUI: Export transparent…
# # Helper: mcf.make_transparent_cutout(...)

# %% [markdown]
# ## Session JSON and export script

# %%
# s.save_state('session.json')
# # later:
# s2 = mcf.McfSession()
# s2.load_state('session.json')
# # mcf.gui(state='session.json')
#
# script = s.export_script()
# open('recreate.py', 'w').write(script)

# %% [markdown]
# Multi-panel **component mosaics** (hero + per-layer strip) are covered in
# [component_mosaic](component_mosaic.py) — that is the preferred multi-panel
# figure path; a mosaic button in the GUI is deferred.

# %%
print('See guide/saving_outputs.md for details.')
