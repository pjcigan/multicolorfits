# Saving outputs

Ways to leave a multicolorfits session with publication-ready products.
Job map: {doc}`../capabilities/index`. Recipes: `mcf.recipes('export')`,
`mcf.recipes('cutout')`, `mcf.recipes('save_rgb')`.

## Rasters

Build a matplotlib figure (GUI **Plot Full Resolution**, or
{func}`~multicolorfits.make_combined_figure`) and save; the **extension
selects the format**. The builder returns a **figure only** — grab
`fig.axes[0]` if you need the axes. Save with `fig.savefig`, not
`plt.savefig` (that often writes a blank canvas).

```python
import multicolorfits as mcf

fig = mcf.make_combined_figure(session)
ax = fig.axes[0]                         # combined WCS axes
fig.savefig('out.png', dpi=300, bbox_inches='tight',
            facecolor=fig.get_facecolor())
fig.savefig('out.pdf', bbox_inches='tight')
```

For image-only stamps (no WCS chrome):

```python
session.apply_bare_plot_style(transparent=True)
fig = mcf.make_combined_figure(session)
fig.savefig('stamp.png', transparent=True, dpi=200, bbox_inches='tight')
```

## RGB FITS

```python
session.save_rgb_fits('combined.fits')           # adds provenance HISTORY
# Low-level / with explicit provenance from a session:
# mcf.save_rgb_fits(path, rgb_array, header)
# mcf.save_combined(path, rgb_array, header, session_info=session)
```

GUI **Save FITS…** can also write annotated combined cubes and per-layer
data/color components with provenance keywords (`MCFTYPE`, `MCFMODE`,
`MCFBLEND`, …).

## Transparent cutouts

For slides with irregular outlines. Session helper (size is the longest
side after crop):

```python
rgba = session.export_transparent_cutout(size=400, savepath='stamp.png')
# GUI: Export transparent…
```

Low-level control of the alpha matte (same kwargs the M74 example uses):

```python
rgba = mcf.make_transparent_cutout(
    combined, crop='auto', pad=0.05,
    alpha_smooth=3, alpha_lo=45, alpha_hi=82, alpha_gamma=0.7)
mcf.save_transparent_cutout(rgba, 'stamp.png')
fig = mcf.preview_cutout_on_backgrounds(rgba)   # checker / dark / light
```

```{image} ../_static/showcase/m74_stamp.png
:alt: M74 transparent-sky cutout stamp
:width: 280px
:align: center
```

Preview the matte on checkerboard / dark / light before shipping:

```{image} ../_static/showcase/m74_cutout_backgrounds_light.png
:class: mcf-plot plot-light
:alt: Cutout previewed over checkerboard, dark, and light backgrounds
```

```{image} ../_static/showcase/m74_cutout_backgrounds_dark.png
:class: mcf-plot plot-dark
:alt: Cutout previewed over checkerboard, dark, and light backgrounds (dark chrome)
```

Worked example with fade shaping (`alpha_gamma`) and deblend helpers:
{doc}`../examples/m74_transparent_cutout`.

## Session JSON and scripts

```python
session.save_state('session.json')
open('recreate.py', 'w').write(session.export_script())
```

`export_script()` replays alignment steps and emits the current PEP 8 names
(not the silent v2 aliases). Tutorial:
{doc}`../tutorials/session_export_script`.

## Multi-panel figures

Prefer {func}`~multicolorfits.figures.make_component_mosaic` for a hero +
per-layer strip (see {doc}`figures_and_mosaics`).  A mosaic export button in
the GUI is deferred; scripting is the supported path.

```python
fig, axes = mcf.make_component_mosaic(session, components='top', ticks='plain')
fig.savefig('mosaic.png', dpi=150, facecolor=fig.get_facecolor())
```

Tutorial: {doc}`../tutorials/save_and_export`.
