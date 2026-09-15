# Capabilities

A job-oriented map of what multicolorfits can do — short starter snippets,
links into the Guide / Examples / API, and the `mcf.recipes('<keyword>')`
query that finds the same idea in the runtime catalog (`llms.txt`).

This is a **thin hub** (not a full feature-gallery codegen lane). Figures
reuse existing showcase / example assets. Deeper prose lives in the Guide;
runnable narratives live under Examples and Tutorials.

```{tip}
Prefer searching by task: ``mcf.recipes('mosaic')``, ``mcf.recipes('beam')``,
``mcf.recipes('gui')``, ``mcf.recipes('export')``.  Full orientation:
``mcf.overview()`` — see {doc}`../api/getting_oriented`.
```

## Coverage matrix

Status of the discoverability surfaces for each capability. Guide method maps,
API “where to look” blurbs, and docstring ``Examples`` on overlays / stack /
WCS / session were deepened after the hub landed; keep this table as the
backlog checklist when adding new features.

| Capability | Guide | Example / tutorial | Recipe keyword | API |
|------------|-------|--------------------|----------------|-----|
| Core colorize → combine | {doc}`../guide/color_compositing` | {doc}`../tutorials/getting_started` | `colorize`, `combine` | {doc}`../api/color_pipeline` |
| One-call / Lab compose | {doc}`../guide/color_compositing` | {doc}`../examples/colorspace_comparison` | `lab`, `combine_layers` | {doc}`../api/compositing` |
| Stretch / Auto levels | {doc}`../guide/intensity_scaling` | — | `suggest`, `levels` | {doc}`../api/color_pipeline` |
| Session editing | {doc}`../guide/session_model` | {doc}`../tutorials/getting_started` | `session`, `set_display` | {doc}`../api/session` |
| Align / north-up / prep | {doc}`../guide/preparing_images` | {doc}`../tutorials/prepare_and_align` | `align`, `prep_layers`, `north` | {doc}`../api/pipeline`, {doc}`../api/wcs_fits` |
| Beam match | {doc}`../guide/preparing_images` | — | `match_beam` | {doc}`../api/wcs_fits` |
| Component mosaic | {doc}`../guide/figures_and_mosaics` | {doc}`../tutorials/component_mosaic` | `mosaic` | {doc}`../api/figures` |
| Overlays | {doc}`../guide/overlays` | {doc}`../examples/ngc602` | `compass`, `scale_bar` | {doc}`../api/overlays` |
| Palettes / preview | {doc}`../guide/palettes` | {doc}`../examples/crab_palette_gallery` | `palette`, `preview_palette` | {doc}`../api/palettes` |
| Transparent cutout | {doc}`../guide/saving_outputs` | {doc}`../examples/m74_transparent_cutout` | `cutout` | {doc}`../api/compositing` |
| Save / export script | {doc}`../guide/saving_outputs` | {doc}`../tutorials/session_export_script` | `export_script`, `save_rgb` | {doc}`../api/session` |
| Browser / Qt / embed GUI | {doc}`../tutorials/gui_web_walkthrough` | {doc}`../tutorials/gui_notebook_embed` | `gui` | {doc}`../api/gui` |

---

## Core pipeline

Stretch → tint → combine (layer-first model).

```{image} ../_static/showcase/ngc602_hero_light.png
:class: mcf-plot plot-light
:alt: NGC 602 combined composite
:width: 420px
```

```{image} ../_static/showcase/ngc602_hero_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 combined composite (dark)
:width: 420px
```

```python
import multicolorfits as mcf

grey = mcf.to_grey_rgb(data, rescalefn='asinh', scaletype='perc', min_max=[1, 99.5])
layer = mcf.colorize_image(grey, '#E4002B', colorintype='hex')
rgb = mcf.combine_multicolor([layer_a, layer_b], gamma=2.2)
```

**See also:** guide {doc}`../guide/color_compositing` · API {py:func}`multicolorfits.to_grey_rgb` · `mcf.recipes('colorize')`

---

## One-call and Lab compositing

```{image} ../_static/showcase/ngc602_rgb_vs_lab_light.png
:class: mcf-plot plot-light
:alt: RGB vs Lab comparison
:width: 480px
```

```{image} ../_static/showcase/ngc602_rgb_vs_lab_dark.png
:class: mcf-plot plot-dark
:alt: RGB vs Lab comparison (dark)
:width: 480px
```

```python
# Arrays → combined RGB in one call
rgb = mcf.combine_layers(
    [data_r, data_g, data_b],
    colors=['#E4002B', '#33CC33', '#0088FF'],
    stretches='asinh', gamma=2.2)

# Same dispatcher the GUI / session use (mode + blend + background)
rgb = mcf.combine_colorized_layers(
    layers, mode='lab', blend='screen', background='black', gamma=2.2)
```

**See also:** {doc}`../guide/color_compositing` · {doc}`../examples/colorspace_comparison` · `mcf.recipes('lab')`

---

## Intensity stretch and Auto levels

```python
rec = mcf.suggest_levels(data)
grey = mcf.to_grey_rgb(
    data, rescalefn=rec['stretch'], scaletype='abs',
    min_max=[rec['vmin'], rec['vmax']])

# Session: known values, or Auto levels on every loaded panel
s.set_display(stretches='asinh', vmins=0.1, vmaxs=12.0)
s.apply_suggested_levels()
```

**See also:** {doc}`../guide/intensity_scaling` · {py:func}`multicolorfits.suggest_levels` · `mcf.recipes('suggest')`

---

## Session editing

```python
s = mcf.McfSession(n_panels=3)
s.load_files(
    [fits_r, fits_g, fits_b],
    colors=['#E4002B', '#33CC33', '#0088FF'],
    labels=['R', 'G', 'B'])
s.compose.combine_mode = 'lab'
s.compose.combine_blend = 'screen'
rgb = s.render_combined()
```

**See also:** {doc}`../guide/session_model` · {py:class}`multicolorfits.McfSession` · `mcf.recipes('session')`

---

## Align, north-up, and prep

```python
# Scripting: list of (data, header)
aligned = mcf.align_stack([(data_r, hdr), (data_g, hdr), (data_b, hdr)],
                          reference=0)
prepared = mcf.prep_layers(
    [(data_r, hdr), (data_g, hdr)],
    north_up=True, oversample=2, crop='overlap')

# Session: same idea on loaded panels
s.align_panels(target='reference', north_up=True, crop='overlap')
```

Needs `pip install "multicolorfits[reproject]"`.

**See also:** {doc}`../guide/preparing_images` · {doc}`../guide/wcs_and_alignment` · {doc}`../tutorials/prepare_and_align` · `mcf.recipes('align')`

---

## Beam matching

```python
# Convolve to a coarser target beam (FWHM arcsec). Explicit *_from when
# the header has no BMAJ/BMIN; Jy/beam maps often want per_beam=True.
smoothed, hdr_out = mcf.match_beam(
    data, header,
    bmaj_to_asec=12.0, bmin_to_asec=12.0, bpa_to_deg=0.0,
    bmaj_from_asec=5.0, bmin_from_asec=5.0, bpa_from_deg=0.0,
    per_beam=False)
```

**See also:** {doc}`../guide/preparing_images` · {py:func}`multicolorfits.match_beam` · `mcf.recipes('match_beam')`

---

## Component mosaic

```{image} ../_static/showcase/ngc602_mosaic_light.png
:class: mcf-plot plot-light
:alt: NGC 602 component mosaic
:width: 480px
```

```{image} ../_static/showcase/ngc602_mosaic_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 component mosaic (dark)
:width: 480px
```

```python
fig, axes = mcf.make_component_mosaic(
    s, components='top', max_per_line=3, ticks='plain')
ax = axes['combined']          # hero RGB — not a session panel
fig.savefig('mosaic.png', dpi=150, facecolor=fig.get_facecolor())
```

**See also:** {doc}`../guide/figures_and_mosaics` · {doc}`../tutorials/component_mosaic` · `mcf.recipes('mosaic')`

---

## Overlays (compass, beam, scale bar)

Canonical guide: {doc}`../guide/overlays` (compose flags vs post-hoc drawing,
which axes is the hero, style kwargs).

```python
# Needs pip install "multicolorfits[overlays]"
if mcf.overlays_available():
    s.compose.show_compass = True
    s.compose.show_scale_bar = True
    s.compose.scale_bar_asec = 30.0
    fig = mcf.make_combined_figure(s)
    # After the fact:
    # mcf.overlays.add_scale_bar(fig.axes[0], s.common_header, length_asec=30)
```

**See also:** {doc}`../api/overlays` · `mcf.recipes('compass')`

---

## Palettes and color-vision preview

```{image} ../_static/showcase/crab_gallery.png
:class: mcf-plot
:alt: Crab palette gallery montage
:width: 520px
```

```python
colors = mcf.suggest_colors(3)                 # CIE LCh auto (GUI Even auto-N)
# colors = mcf.colors_for_hue_pattern('triad', n=3)
# colors = mcf.colors_from_hsv(3)              # classical HSV wheel
res = mcf.preview_palette(colors, mode='lab', blend='screen',
                          background='black')
# plt.show() or res.show() — figure is pyplot-managed
res.fig.savefig('palette_preview.png', dpi=120,
                facecolor=res.fig.get_facecolor())
# From a session: s.preview_palette() / s.apply_palette('pob')
```

**See also:** {doc}`../guide/palettes` · {doc}`../examples/crab_palette_gallery` · {doc}`../api/palettes` · `mcf.recipes('palette')`

---

## Transparent cutouts

```{image} ../_static/showcase/m74_stamp.png
:class: mcf-plot
:alt: M74 transparent stamp
:width: 280px
```

```python
rgba = mcf.make_transparent_cutout(
    combined, crop='auto', pad=0.05,
    alpha_smooth=3, alpha_lo=45, alpha_hi=82, alpha_gamma=0.7)
mcf.save_transparent_cutout(rgba, 'stamp.png')
# Session helper: s.export_transparent_cutout(size=512, savepath='stamp.png')
```

**See also:** {doc}`../examples/m74_transparent_cutout` · {doc}`../guide/saving_outputs` · `mcf.recipes('cutout')`

---

## Save outputs and export script

```python
s.save_rgb_fits('combined_rgb.fits')
script = s.export_script()          # standalone recreate script
open('recreate.py', 'w').write(script)

# Provenance-annotated FITS from a finished array:
mcf.save_combined('combined_rgb.fits', combined, header, session_info=s)
```

**See also:** {doc}`../guide/saving_outputs` · {doc}`../tutorials/session_export_script` · `mcf.recipes('export')`

---

## Browser, Qt, and notebook GUIs

```{image} ../_static/gui/gui_web_walkthrough_light.png
:class: mcf-plot plot-light
:alt: Browser GUI walkthrough
:width: 480px
```

```{image} ../_static/gui/gui_web_walkthrough_dark.png
:class: mcf-plot plot-dark
:alt: Browser GUI walkthrough (dark)
:width: 480px
```

```python
# After building or loading a session (do not call these in batch scripts):
# mcf.gui(session=s)                 # browser — pip install multicolorfits[web]
# mcf.gui_embed(session=s, height=900)
# mcf.gui_qt(session=s)              # desktop — pip install multicolorfits[qt]
# mcf.gui(files=[fits_r, fits_g, fits_b], colors=[...])
```

**See also:** {doc}`../tutorials/gui_web_walkthrough` · {doc}`../tutorials/gui_notebook_embed` · {doc}`../api/gui` · `mcf.recipes('gui')`
