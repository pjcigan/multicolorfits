# Saving outputs

Ways to leave a multicolorfits session with publication-ready products.

## Rasters

Build a matplotlib figure (GUI **Plot Full Resolution**, or
`make_combined_figure`) and save; the **extension selects the format**:

```python
from multicolorfits.figures import make_combined_figure

fig, ax = make_combined_figure(session)
fig.savefig('out.png', dpi=300, bbox_inches='tight')
fig.savefig('out.pdf', bbox_inches='tight')
```

For image-only stamps (no WCS chrome):

```python
session.apply_bare_plot_style(transparent=True)
fig, ax = make_combined_figure(session)
fig.savefig('stamp.png', transparent=True, dpi=200, bbox_inches='tight')
```

## RGB FITS

```python
session.save_rgb_fits('combined.fits')           # adds provenance HISTORY
# mcf.save_rgb_fits(path, rgb_array, header)       # low-level
```

GUI **Save FITS…** can also write annotated combined cubes and per-layer
data/color components with provenance keywords (`MCFTYPE`, `MCFMODE`,
`MCFBLEND`, …).

## Transparent cutouts

For slides with irregular outlines:

```python
rgba = session.export_transparent_cutout(size=400)
# GUI: Export transparent…
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
(not the silent v2 aliases).

## Multi-panel figures

Prefer {func}`~multicolorfits.figures.make_component_mosaic` for a hero +
per-layer strip (see {doc}`figures_and_mosaics`).  A mosaic export button in
the GUI is deferred; scripting is the supported path.

Tutorial: {doc}`../tutorials/save_and_export`.
