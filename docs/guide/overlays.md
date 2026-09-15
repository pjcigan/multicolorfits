# Overlays and annotations

Scale bars, compass roses, beams, color swatches, and band labels are drawn
on the **combined matplotlib axes**. They are not extra entries in
`session.panels`. That list is only the input layers.

Job map: {doc}`../capabilities/index`. API: {doc}`../api/overlays`.
Recipe: `mcf.recipes('compass')` / `mcf.recipes('scale_bar')`.

Two ways to add them:

1. Set flags on `session.compose`, then build the figure. The builders read
   those flags.
2. Build the figure first, grab the combined axes, and draw onto it. Use
   this when the figure already exists, or when a mosaic was built with
   `ticks='plain'` (that mode skips the channel legend).

Compass, beam, and scale bar need the optional extra:

```bash
pip install "multicolorfits[overlays]"
```

Check with `mcf.overlays_available()`. Legend, combo swatch, and band labels
are in core and do not need that extra. API reference:
{doc}`../api/overlays`.

---

## Which axes is the combined image?

```python
import multicolorfits as mcf

fig, axes = mcf.make_component_mosaic(s, components='top', ticks='plain')
ax = axes['combined']            # hero: the composited RGB image
strip = axes['components']       # one axes per loaded layer, same order
hdr = s.common_header            # WCS / pixel scale for overlays
```

`len(s.panels)` does not grow when you make the mosaic. There is no
"combined panel" object.

A single combined figure does not return the axes:

```python
fig = mcf.make_combined_figure(s)
ax = fig.axes[0]                 # grab this before adding corner insets
```

`setup_combined_axes` returns the axes if you would rather not dig it out
of the figure:

```python
from matplotlib.figure import Figure

fig = Figure(figsize=(7, 7), facecolor='white')
ax = mcf.setup_combined_axes(fig, s)
```

Either figure is a standalone `Figure`, not the pyplot current figure. Save
with `fig.savefig(...)`, not `plt.savefig` (that writes a blank canvas).

`overlays='hero'` (the mosaic default) draws compass / beam / scale bar on
the hero only, and only if the matching `compose` flags are already on.
`overlays='none'` leaves the hero clean so you can add them yourself.

---

## Set them before the figure is built

```python
s.compose.show_legend = True
s.compose.legend_loc = 'upper right'
s.compose.show_combo_swatch = True
s.compose.combo_swatch_labels = True
s.compose.combo_swatch_loc = 'lower left'
s.compose.show_band_labels = True
s.compose.band_labels_loc = 'upper left'

s.compose.show_compass = True
s.compose.compass_loc = 'upper right'     # corner name
s.compose.show_scale_bar = True
s.compose.scale_bar_asec = 1.0            # arcsec; 0 = auto (~15% of the field)
s.compose.scale_bar_loc = 2               # integer, not a corner name — see below
s.compose.scale_bar_color = 'white'
s.compose.show_beam = True                # needs BMAJ / BMIN in the header
s.compose.beam_style = 'ellipse'          # or 'crosshair', 'crosshairgrid'

fig, axes = mcf.make_component_mosaic(s, ticks='minimal', overlays='hero')
fig.savefig('mosaic.jpg', dpi=200, bbox_inches='tight',
            facecolor=fig.get_facecolor())
```

`ticks='plain'` hides ticks and also skips the channel legend, even when
`show_legend` is True. The combo swatch is still drawn. Compass, beam, and
scale bar follow `overlays=` (`'hero'` or `'none'`), not the tick mode.

If skyplothelper is missing, the compass / beam / scale-bar flags warn and
are skipped. The rest of the figure still builds.

---

## Add them after the figure exists

This is the path when you already called `make_component_mosaic` and want
to decorate the hero, including a legend that `ticks='plain'` omitted.

```python
from multicolorfits.figures import add_band_labels, add_swatch_inset

fig, axes = mcf.make_component_mosaic(
    s, components='top', ticks='plain', overlays='none')
ax = axes['combined']
hdr = s.common_header

add_swatch_inset(
    ax, s.render_combo_swatch(),
    loc='lower left', show_labels=True, inset_scale=0.22)
add_band_labels(ax, s.legend_entries(), loc='upper left')

if mcf.overlays_available():
    mcf.overlays.add_compass(ax, loc='upper right', color='white')
    mcf.overlays.add_scale_bar(
        ax, hdr, length_asec=1.0, loc=2, color='white')
    # Optical / JWST mosaics often have no beam cards; this warns and skips.
    mcf.overlays.add_beam(ax, hdr, loc='lower left', style='ellipse')

fig.savefig('mosaic.jpg', dpi=200, bbox_inches='tight',
            facecolor=fig.get_facecolor())
```

One call for the skyplothelper trio, using the same arguments as the
compose flags:

```python
mcf.overlays.apply_overlays(
    ax, hdr,
    compass=True, compass_loc='upper right',
    scale_bar=True, scale_bar_asec=1.0, scale_bar_loc=2,
    beam=False,
)
```

`apply_compose_overlays(ax, hdr, s.compose)` reads the flags already stored
on the session and draws those. Useful if you built the mosaic with
`overlays='none'` and then flipped the flags.

The strip axes (`axes['components']`) are ordinary WCS axes too, so the
same helpers work there. A scale bar on every component is usually noise;
the hero is the one readers measure.

---

## What each annotation is

### Channel legend and band labels

A legend is a matplotlib key: color patch plus the panel `label` (or
`Image N` if the label is empty). Band labels are the same names drawn in
the layer color, chained in a corner, with a thin contrasting outline so
they stay readable on the image.

```python
s.legend_entries()    # [(hex, label), ...] for loaded panels
```

### Combo swatch

Overlapping circles, one per layer, mixed with the **same** combine mode,
blend, gamma, and background as the hero. Overlap hues in the swatch match
the image. Lab / RYB / CMYK overlaps are not a generic RGB diagram.

`s.render_combo_swatch()` returns the image and circle positions.
`add_swatch_inset` places that image in a corner with no frame.
`inset_scale` is the inset size as a fraction of the axes (default 0.24).
`label_offset` pushes each name into that circle's outer lobe (default 0.62).

### Compass

North and east arrows that follow the WCS, not the pixel axes. A north-up
image still gets a compass; a rotated field gets arrows that are not aligned
with the plot border.

`loc` is a corner name (`'upper right'`, `'lower left'`, …) or an
`(x, y)` pair in axes fractions. `length` is the arrow length as a fraction
of the axes (default 0.08). A named corner places the tail at
`pad + length`, so the arrow length is most of the inset. `pad` is extra
inset beyond that, in axes fraction. The multicolorfits default is 0.02
(skyplothelper's 0.05 sits the letters near a quarter of the way into the
frame). `label_offset` (default 1.3, in arrow lengths) is how far the
letters sit past the heads. `fontsize` is in points. `color` is the arrow
and letter ink. `stroke_color` / `stroke_lw` outline them; see
{ref}`overlay-stroke`. `north_label` / `east_label` rename the letters.

The axes must be a WCS axes. Both figure builders project the common header
when it has a WCS, so `axes['combined']` qualifies.

### Scale bar

An angular bar from the header pixel scale (`CDELT` / `CD` / `PC`). Pass
`length_asec` in arcseconds. Omit it, or pass 0, and the length is about
15% of the field width, rounded to 1 / 2 / 5 × 10ⁿ. The label is `30"` or
`2'` from that length unless you pass `label=`.

**Corner codes are not the same as the compass.** The scale bar uses
matplotlib's location integers:

| `scale_bar_loc` | Corner |
|-----------------|--------|
| 1 | upper right |
| 2 | upper left |
| 3 | lower left |
| 4 | lower right (default) |

`color` is the bar; `stroke_color` and `stroke_lw` outline it so a white
bar stays visible on a bright background (defaults `'k'` and 1.75 — see
{ref}`overlay-stroke`). Compose fields: `scale_bar_color` (default
`'white'`), `scale_bar_stroke_color` (`'black'`), `scale_bar_stroke_lw`
(1.75).

The header has to carry a pixel scale. `s.common_header` is the first
loaded panel's header; that is the right one after the layers share a grid.

### Beam

The synthesized beam from `BMAJ`, `BMIN`, and `BPA`. Radio cubes usually
have these. HST / JWST science headers often do not. Without the cards,
`add_beam` warns and draws nothing — the figure is otherwise unchanged.

`style` is `'crosshair'` (default), `'ellipse'`, or `'crosshairgrid'`.
`ec` / `lw` / `fc` style the ellipse. `stroke_color` is off until set
({ref}`overlay-stroke`). `loc` is a corner name, as for the compass.
`anchored=False` draws the ellipse at the beam centre in data coordinates
instead of pinning it to a corner.

---

## Without a session

`plot_combined_rgb` takes the same overlay switches and saves the figure
itself. It does not return the axes, so it is the "set them and save" path,
not the "add them afterwards" path:

```python
mcf.plot_combined_rgb(
    rgb, header, 'NGC 602', 'out.png',
    show_compass=True, compass_loc='upper right',
    show_scale_bar=True, scale_bar_asec=30, scale_bar_loc=4,
)
```

---

(overlay-stroke)=
## Color and stroke

`stroke_color` / `stroke_lw` already pass through to skyplothelper. You do
not need a separate mcf stroke helper, and you do not have to set matplotlib
path effects unless you want to change the stroke after the artist exists.

```python
mcf.overlays.add_compass(
    ax, loc='upper right', color='white',
    stroke_color='black', stroke_lw=2)
mcf.overlays.add_scale_bar(
    ax, hdr, length_asec=1.0, loc=2,
    color='white', stroke_color='black', stroke_lw=2)
mcf.overlays.add_beam(
    ax, hdr, loc='lower left', style='ellipse',
    ec='white', stroke_color='black')
```

Scale-bar stroke defaults to a thin black outline (`stroke_color='k'`,
`stroke_lw=1.75`). Compass stroke defaults to the axes facecolor. Beam
stroke is off until you set `stroke_color`. `stroke_color=None` on the
scale bar turns that default outline off. The scale-bar label has no
`fontsize` keyword — pass `fontproperties` (or `prop`):

There is no `fontsize` keyword. A dict works:

```python
mcf.overlays.add_scale_bar(
    ax, hdr, length_asec=1.0, loc=2,
    fontproperties={'size': 8})
```

Through the batch helper, styling lives in a per-overlay dict so `color=`
is not ambiguous:

```python
mcf.overlays.apply_overlays(
    ax, hdr, compass=True, scale_bar=True,
    compass_kwargs=dict(color='white', stroke_color='black'),
    scale_bar_kwargs=dict(color='white', stroke_lw=2),
)
```

The session path only stores scale-bar stroke
(`s.compose.scale_bar_color`, `scale_bar_stroke_color`,
`scale_bar_stroke_lw`). Compass and beam stroke are kwargs on the
`add_*` calls, or matplotlib `path_effects` on the returned artists.

Other forwarded keywords are listed on
{func}`~multicolorfits.overlays.add_compass`,
{func}`~multicolorfits.overlays.add_scale_bar`, and
{func}`~multicolorfits.overlays.add_beam` (size, font, corner padding,
arrow-head size, beam edge/face, crosshair style). Anything not named there
is an escape hatch into the underlying skyplothelper artist, not a second
mcf API.

---

## Further skyplothelper

`mcf.overlays.skyplothelper()` returns the skyplothelper module (and raises
`ImportError` with the install line if it is missing). Offset-coordinate
frames and other specialized axes stay there:

```python
sph = mcf.overlays.skyplothelper()
wcs = mcf.overlays.offset_coord_wcs(hdr, center)
sph.add_coord_overlay(ax, wcs)
```

`mcf.overlays.make_offset_figure(center)` is the one-call offset-field
figure. Prefer that module for graticules and globe frames rather than
reimplementing them here.
