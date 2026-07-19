# Figures and mosaics

## Combined WCS figure

```python
from multicolorfits.figures import make_combined_figure

rgb = session.render_combined()
fig, ax = make_combined_figure(session, combined=rgb)
```

Compose flags on `session.compose` control ticks, title, legend, combo swatch,
band labels, and optional skyplothelper overlays.  Legend / swatch concepts:
{doc}`session_model`.

```{image} ../_static/showcase/ngc602_hero_light.png
:class: mcf-plot plot-light
:alt: NGC 602 Lab/screen combined figure with WCS chrome
```

```{image} ../_static/showcase/ngc602_hero_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 Lab/screen combined figure (dark chrome)
```

## Component mosaic

Hero combined panel plus a strip of colorized components on a **shared WCS**:

```python
import multicolorfits as mcf

fig, axes = mcf.make_component_mosaic(
    session,
    components='top',     # or bottom / left / right
    max_per_line=3,
    ticks='plain',        # plain | minimal | full
    overlays='hero',
)
# axes['combined'], axes['components']
```

Use this when you want a figure that explains each layer, not only the blend.

### Layout examples

**Components on top**, wrapping at three per row (`max_per_line=3`):

```{image} ../_static/mosaics/mosaic_top_max3_light.png
:class: mcf-plot plot-light
:alt: Component mosaic with strip on top, max 3 per row
```

```{image} ../_static/mosaics/mosaic_top_max3_dark.png
:class: mcf-plot plot-dark
:alt: Component mosaic with strip on top (dark chrome)
```

**Tighter strip** (`max_per_line=2`):

```{image} ../_static/mosaics/mosaic_top_max2_light.png
:class: mcf-plot plot-light
:alt: Component mosaic, max 2 per row
```

```{image} ../_static/mosaics/mosaic_top_max2_dark.png
:class: mcf-plot plot-dark
:alt: Component mosaic, max 2 per row (dark chrome)
```

**Components on the left** of the hero:

```{image} ../_static/mosaics/mosaic_left_max3_light.png
:class: mcf-plot plot-light
:alt: Component mosaic with strip on the left
```

```{image} ../_static/mosaics/mosaic_left_max3_dark.png
:class: mcf-plot plot-dark
:alt: Component mosaic with strip on the left (dark chrome)
```

Tutorial notebook: {doc}`../tutorials/component_mosaic`.
Local regen scratch → `_static/mosaics/` via `python docs/make_docs_figures.py`.

## Overlays (compass, beam, scale bar)

Optional chrome from skyplothelper (`pip install "multicolorfits[overlays]"`):

```python
s.compose.show_compass = True
s.compose.show_beam = True          # needs BMAJ/BMIN (etc.) in the header
s.compose.show_scale_bar = True
# locations: compose.compass_loc, beam_loc, scale_bar_loc, …

fig, ax = make_combined_figure(s)
# or on a mosaic hero only:
fig, axes = mcf.make_component_mosaic(s, overlays='hero')
```

Check availability with `mcf.overlays_available()`.  Low-level helpers:
{doc}`../api/overlays`.  Export Script emits the same toggles when they are on.
