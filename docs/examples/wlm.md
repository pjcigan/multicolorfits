# WLM

WLM (Wolf–Lundmark–Melotte) is a nearby low-metallicity dwarf. Download HI /
optical / UV products from the
[LITTLE THINGS NRAO page](https://science.nrao.edu/science/surveys/littlethings/data/wlm.html).

This walkthrough uses **HI**, **V-band**, and **near-UV**. The goal is a
balanced three-color image that roughly matches the NRAO thumbnail aesthetic
(which was finished in Photoshop on a simpler RGB base).

Companion notebook: {doc}`../tutorials/wlm_align_reproject`. Concepts:
{doc}`../guide/color_compositing`, {doc}`../guide/session_model`.

```{note}
Layer colors here are red / yellow / blue — that is **not** the same as the
subtractive ``combine_mode='ryb'`` paint-mix compositor.
```

---

## Recommended — crop, align, session

```python
import multicolorfits as mcf

# Load originals; >2D cubes are squeezed automatically on PanelState.load_fits.
s = mcf.McfSession(n_panels=3)
s.load_files(
    ['./WLM_NA_X0_P_R.FITS', './wlmv.fits', './wlmncut.fit'],
    colors=['#994242', '#FFF9DB', '#1773E9'],
    labels=['HI', 'V', 'NUV'],
)

# Prefer align_stack / GUI Align when grids differ (needs [reproject]).
# After an interactive crop in the GUI, Save Session and reload — or crop in
# script first (classic path below), then:
# s.align_panels(target='reference')

s.compose.combine_mode = 'lab'
s.compose.combine_blend = 'screen'
s.compose.gamma = 2.2
s.compose.title = 'WLM — HI, V, NUV'

# mcf.gui(session=s)
rgb = s.render_combined()
s.save_rgb_fits('./wlm_lab.fits')
print(s.export_script()[:500], '...')
```

---

## Classic scripting pipeline

### Load and tidy HI header

```python
import numpy as np
import astropy.io.fits as pyfits
import multicolorfits as mcf

wlm_hidat, wlm_hihdr = pyfits.getdata('./WLM_NA_X0_P_R.FITS', header=True)
wlm_vdat, wlm_vhdr = pyfits.getdata('./wlmv.fits', header=True)
wlm_nuvdat, wlm_nuvhdr = pyfits.getdata('./wlmncut.fit', header=True)

# Explicit squeeze is still useful when working with raw arrays:
wlm_hidat = np.squeeze(wlm_hidat)
wlm_hihdr = mcf.force_header_2d(wlm_hihdr)
wlm_hihdr['RADESYS'] = 'FK5'   # often missing on this HI product

wlm_hihdr_simple = mcf.make_simple_header(wlm_hihdr, radesys='FK5')
```

### Crop HI as master grid

```python
cropcenter_coords = mcf.pixel_to_sky(wlm_hihdr_simple, 512, 512)
cropwidth_asec = 600.  # half-width from center

wlm_hicropdat, wlm_hicrophdr = mcf.crop_image_sky(
    wlm_hidat, wlm_hihdr_simple, cropcenter_coords, cropwidth_asec,
    savenew='./wlm_hicrop.fits', overwrite=True)
```

### Reproject companions (or `align_stack`)

```python
# Explicit per-layer reproject (v2-style):
wlm_vcropdat = mcf.reproject_image(
    wlm_vdat, mcf.make_simple_header(wlm_vhdr), wlm_hicrophdr)
pyfits.writeto('./wlm_vcrop.fits', wlm_vcropdat, wlm_hicrophdr, overwrite=True)

wlm_nuvcropdat = mcf.reproject_image(
    wlm_nuvdat, mcf.make_simple_header(wlm_nuvhdr), wlm_hicrophdr)
pyfits.writeto('./wlm_nuvcrop.fits', wlm_nuvcropdat, wlm_hicrophdr, overwrite=True)

# Preferable when you already have arrays + headers in a list:
# aligned, hdr = mcf.align_stack(
#     [wlm_hidat, wlm_vdat, wlm_nuvdat],
#     [wlm_hihdr_simple, wlm_vhdr, wlm_nuvhdr],
#     reference=0)
```

Reload the cropped products for GUI use:

```python
# mcf.gui(files=['./wlm_hicrop.fits', './wlm_vcrop.fits', './wlm_nuvcrop.fits'],
#         colors=['#994242', '#FFF9DB', '#1773E9'],
#         labels=['HI', 'V', 'NUV'])
```

---

## Colorize and combine

Leave `checkscale=False` (default) for non-interactive runs; use the GUI or
percentile / absolute limits for level tweaking.

```python
hi_greyRGB = mcf.to_grey_rgb(
    wlm_hicropdat, rescalefn='asinh', scaletype='perc',
    min_max=[1., 99.9], gamma=2.2)

v_greyRGB = mcf.to_grey_rgb(
    wlm_vcropdat, rescalefn='asinh', scaletype='abs',
    min_max=[3650., 4800.], gamma=2.2)

nuv_greyRGB = mcf.to_grey_rgb(
    wlm_nuvcropdat, rescalefn='asinh', scaletype='perc',
    min_max=[20, 99.9], gamma=2.2)

hi_red = mcf.colorize_image(hi_greyRGB, '#994242', colorintype='hex', gammacorr_color=2.2)
v_yellow = mcf.colorize_image(v_greyRGB, '#FFF9DB', colorintype='hex', gammacorr_color=2.2)
nuv_blue = mcf.colorize_image(nuv_greyRGB, '#1773E9', colorintype='hex', gammacorr_color=2.2)

# Classic RGB sum, or Lab/screen for softer overlaps:
wlm_custom = mcf.combine_multicolor([hi_red, v_yellow, nuv_blue], gamma=2.2)
# wlm_lab = mcf.combine_multicolor_colorspace(
#     [hi_red, v_yellow, nuv_blue], colorspace='lab', blend='screen', gamma=2.2)

mcf.plot_combined_rgb(
    wlm_custom, wlm_hicrophdr, 'WLM — HI (red), V (yellow), NUV (blue)',
    './WLM_testplot.jpg', tickcolor='0.6', labelcolor='k', facecolor='w',
    minorticks=True, dpi=150)
```

```{image} ../_static/examples/WLM_testplot.jpg
:class: mcf-plot
:alt: WLM HI, V, and NUV custom multicolor
```

### Pure RGB comparison

```python
wlm_pureRGB = np.dstack([
    hi_greyRGB[:, :, 0], v_greyRGB[:, :, 0], nuv_greyRGB[:, :, 0]])

mcf.compare_multicolor_vs_rgb(
    wlm_pureRGB, wlm_custom, wlm_hicrophdr,
    'Custom Multicolor: red/yellow/blue', 'WLM', './wlm_compare.jpg',
    tickcolor='0.6', supy=.75)
```

```{image} ../_static/examples/wlm_compare.jpg
:class: mcf-plot
:alt: Simple RGB stack versus custom red/yellow/blue colors on WLM
```

### Save RGB FITS

```python
mcf.save_rgb_fits('./wlm_custom.fits', wlm_custom, wlm_hicrophdr)
# Session equivalent (adds provenance HISTORY):
# s.save_rgb_fits('./wlm_lab.fits')
```
