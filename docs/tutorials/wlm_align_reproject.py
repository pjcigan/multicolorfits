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
# # WLM — crop, align, combine
#
# Data: [LITTLE THINGS / WLM](https://science.nrao.edu/science/surveys/littlethings/data/wlm.html)
# (HI, V, NUV).  Full narrative: [examples/wlm](../examples/wlm.md).
#
# Layer colors red / yellow / blue are **not** the same thing as subtractive
# `combine_mode='ryb'` — see the [color compositing guide](../guide/color_compositing.md).

# %%
import numpy as np
import astropy.io.fits as pyfits
import multicolorfits as mcf

# %% [markdown]
# ## Squeeze and crop HI as the master grid
#
# `PanelState.load_fits` / `McfSession.load_files` squeeze phantom axes
# automatically.  Explicit squeeze is only needed when you load arrays by hand.

# %%
# hi, hihdr = pyfits.getdata('./WLM_NA_X0_P_R.FITS', header=True)
# v, vhdr = pyfits.getdata('./wlmv.fits', header=True)
# nuv, nuvhdr = pyfits.getdata('./wlmncut.fit', header=True)
#
# hi = np.squeeze(hi)
# hihdr = mcf.force_header_2d(hihdr)
# hihdr['RADESYS'] = 'FK5'
# hi_simple = mcf.make_simple_header(hihdr, radesys='FK5')
#
# center = mcf.pixel_to_sky(hi_simple, 512, 512)
# hi_c, hi_chdr = mcf.crop_image_sky(
#     hi, hi_simple, center, 600., savenew='./wlm_hicrop.fits', overwrite=True)

# %% [markdown]
# ## Align companions
#
# Prefer `align_stack` (or GUI **Align layers…**) when `[reproject]` is installed.
# The classic `reproject_image` calls remain available.

# %%
# # Option A — stack helper
# datas, hdr = mcf.align_stack(
#     [hi_c, v, nuv],
#     [hi_chdr, mcf.make_simple_header(vhdr), mcf.make_simple_header(nuvhdr)],
#     reference=0)
#
# # Option B — explicit reproject onto cropped HI
# # v_c = mcf.reproject_image(v, mcf.make_simple_header(vhdr), hi_chdr)
# # nuv_c = mcf.reproject_image(nuv, mcf.make_simple_header(nuvhdr), hi_chdr)
# # pyfits.writeto('./wlm_vcrop.fits', v_c, hi_chdr, overwrite=True)
# # pyfits.writeto('./wlm_nuvcrop.fits', nuv_c, hi_chdr, overwrite=True)

# %% [markdown]
# ## Session compose (Lab recommended)

# %%
# s = mcf.McfSession(n_panels=3)
# s.load_files(
#     ['./wlm_hicrop.fits', './wlm_vcrop.fits', './wlm_nuvcrop.fits'],
#     colors=['#994242', '#FFF9DB', '#1773E9'],
#     labels=['HI', 'V', 'NUV'],
# )
# for p, stretch, lo, hi in [
#     (s.panels[0], 'asinh', None, None),
# ]:
#     p.stretch = 'asinh'
# # Or tweak levels in the GUI:
# # mcf.gui(session=s)
#
# s.compose.combine_mode = 'lab'
# s.compose.combine_blend = 'screen'
# s.compose.gamma = 2.2
# s.compose.title = 'WLM — HI, V, NUV'
# s.save_rgb_fits('./wlm_lab.fits')
# s.save_state('./wlm_session.json')
# print(s.export_script()[:400], '...')

# %% [markdown]
# ## Classic colorize path (legacy RGB sum)
#
# Use absolute/percentile limits via `to_grey_rgb`.  Keep
# `checkscale=False` for non-interactive runs.

# %%
# hi_g = mcf.to_grey_rgb(hi_c, rescalefn='asinh', scaletype='perc',
#                             min_max=[1., 99.9], gamma=2.2)
# # ... colorize + combine_multicolor / combine_multicolor_colorspace ...

# %%
print('Uncomment after downloading LITTLE THINGS WLM products.')
