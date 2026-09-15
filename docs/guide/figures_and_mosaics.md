# Figures and mosaics

Job snippets: {doc}`../capabilities/index`. Recipes: `mcf.recipes('mosaic')`,
`mcf.recipes('figure')`.

## Combined WCS figure

`session.panels` are the input layers only. The combined image is an RGB
array (`session.render_combined()`) drawn onto a matplotlib axes when you
build a figure. That axes is what you annotate.

```python
fig = mcf.make_combined_figure(session)
ax = fig.axes[0]                 # combined axes; grab it before adding insets
# ax = mcf.setup_combined_axes(fig, session)  # same axes, returned directly
fig.savefig('combined.png', dpi=200, bbox_inches='tight',
            facecolor=fig.get_facecolor())
```

`make_combined_figure` returns the figure only, and that figure is not the
pyplot current figure, so `plt.savefig` writes a blank canvas.

Compose flags on `session.compose` (legend, combo swatch, band labels,
compass, beam, scale bar) are applied while the figure is built. Set them
first, or draw onto `ax` afterwards — see {doc}`overlays`.

```{image} ../_static/showcase/ngc602_hero_light.png
:class: mcf-plot plot-light
:alt: NGC 602 Lab/screen combined figure with WCS chrome
```

```{image} ../_static/showcase/ngc602_hero_dark.png
:class: mcf-plot plot-dark
:alt: NGC 602 Lab/screen combined figure (dark chrome)
```

## Component mosaic

Hero combined panel plus a strip of colorized components on a **shared WCS**.
Each strip panel is that layer as it was prepared for the composite (its
stretch, limits, color, and display gamma), not a reload of the raw FITS:

```python
import multicolorfits as mcf

fig, axes = mcf.make_component_mosaic(
    session,
    components='top',     # or bottom / left / right
    max_per_line=3,
    ticks='plain',        # plain | minimal | full
    overlays='hero',      # compass / beam / scale bar on the hero only
)
ax = axes['combined']         # composited image — not s.panels[i]
# axes['components'][i]       # strip panel i, same order as the loaded layers
fig.savefig('mosaic.png', dpi=200, bbox_inches='tight',
            facecolor=fig.get_facecolor())
```

`s.panels` stays length 3 (or however many files you loaded). Nothing in
that list is the blended result. Overlays and a color swatch go on
`axes['combined']`. With `ticks='plain'` the channel legend is skipped
even if `show_legend` is set; add it after the fact, or use
`ticks='minimal'`. Full annotation examples: {doc}`overlays`.

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
Regen layout figures (dual light/dark, no chrome remap)::

    python docs/make_docs_figures.py --mosaics-only

## Overlays and other annotations

Compass, beam, scale bar, color swatch, and band labels are drawn on the
combined axes, not stored as another panel. Set the `session.compose` flags
before building the figure, or add them afterwards to `axes['combined']` /
`fig.axes[0]`.

{doc}`overlays` covers both paths, including the `[overlays]` extra
(skyplothelper) and how scale-bar corners differ from compass corners.
