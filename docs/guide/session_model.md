# Session model

Both GUIs are thin views over a shared in-memory controller:
{class}`~multicolorfits.session.McfSession`.  Scripts can use the same object
without opening a window. Job-oriented snippets:
{doc}`../capabilities/index`.

## Pieces

| Object | Role |
|--------|------|
| `PanelState` | One image layer: data, header/`McfHeader`, stretch, limits, color, label, smooth |
| `ComposeState` | Combined-figure settings: mode, blend, background, gamma, ticks, legend, overlays, … |
| `McfSession` | `panels[]` + `compose`; render, align, save, export |

## Method map (I want X → call Y)

| Task | Method / attribute |
|------|--------------------|
| Load FITS (optional colors, stretches, limits) | `load_files(...)` |
| In-memory arrays | `panels[i].set_data(data, header)` |
| Known stretch / vmin / vmax | `set_display(...)` |
| Auto levels (histogram start) | `apply_suggested_levels()` |
| Lab / RGB / blend / background | `compose.combine_mode`, `combine_blend`, `combine_background` |
| Align onto one grid | `align_panels(...)` — needs `[reproject]` |
| Combined RGB array | `render_combined()` / `render_combined(preview=True)` |
| Hero WCS figure | `make_combined_figure(s)` → `fig.axes[0]` |
| Hero + component strip | `make_component_mosaic(s)` or `s.plot_component_mosaic(...)` |
| Compass / scale bar / beam flags | `compose.show_*` — see {doc}`overlays` |
| Palette / CVD check / preview | `apply_palette`, `colorblind_report`, `preview_palette` |
| Transparent stamp | `export_transparent_cutout(...)` |
| RGB FITS / session JSON / recreate script | `save_rgb_fits`, `save_state` / `load_state`, `export_script` |
| Open a GUI on this session | `mcf.gui(session=s)`, `gui_embed`, `gui_qt` |

API signatures: {doc}`../api/session`. Recipes: `mcf.recipes('session')`,
`mcf.recipes('export')`, `mcf.recipes('gui')`.

```python
import multicolorfits as mcf

s = mcf.McfSession()                 # default 4 empty panels
s = mcf.McfSession(n_panels=3)

s.load_files(['a.fits', 'b.fits'], colors=['#f00', '#0ff'], labels=['A', 'B'])
s.apply_suggested_levels()          # Auto levels on every loaded panel
# or set values you already chose (a scalar broadcasts):
s.set_display(stretches='asinh', vmins=[0.2, 0.1], vmaxs=[12, 8])
# or characterize in-memory layers first (quiet; returns colors + levels):
# report = mcf.describe_images({'A': (data_a, hdr_a), 'B': (data_b, hdr_b)})
# s.load_files(paths, colors=report['colors'], labels=report['names'],
#              stretches=report['stretches'],
#              vmins=report['vmins'], vmaxs=report['vmaxs'])
s.panels[0].stretch = 'log'         # one panel
s.panels[0].set_limits(0.5, 20)
s.compose.combine_mode = 'lab'
s.compose.combine_blend = 'screen'

rgb = s.render_combined()
```

Stretch and limit suggestions: {doc}`intensity_scaling` (`describe_image`,
`describe_images`, Auto levels).

## Dynamic panel count

Defaults match the GUI (**4** tabs).  Grow or shrink in scripts or via
**+ Add panel** / **Remove** in either GUI (bounds: 1–16):

```python
from multicolorfits import MAX_PANELS, MIN_PANELS  # 16, 1

idx = s.add_panel(color='#88FF88', label='extra')
s.remove_panel(idx)
s.set_n_panels(6)
```

`load_files` grows the list when given more paths than current slots.
`load_state` restores the saved slot count (clamped to `MAX_PANELS`).

## Figure chrome — legend and combo swatch

Publication figures usually need to say which color is which band, and whether
overlap hues are honest for the current combine mode:

```python
s.compose.show_legend = True          # channel → color key
s.compose.show_combo_swatch = True    # overlapping circles using the same mixer
s.compose.show_band_labels = True     # optional text on the image
```

The combo swatch calls the same `combine_colorized_layers` path as the main
image, so Lab / RYB / RGB overlaps in the legend match the figure. These
flags are read when a figure is built. They are not a fourth panel —
`s.panels` stays the input layers. To draw them onto an existing axes, see
{doc}`overlays`.

```{image} ../_static/compositing/swatch_legend_demo_light.png
:class: mcf-plot plot-light
:alt: Channel legend and combo swatch on a combined figure
```

```{image} ../_static/compositing/swatch_legend_demo_dark.png
:class: mcf-plot plot-dark
:alt: Channel legend and combo swatch (dark chrome)
```

Real-data example with legend + swatch: {doc}`../examples/ngc602`.

## Fast preview vs full resolution

Interactive edits on large FITS should not recompute a full float64 WCS figure
on every slider move:

```python
preview = s.render_combined(preview=True)   # float32, longest side ≤ 1024 by default
final = s.render_combined()                 # export / Plot Full Resolution path
```

Both GUIs keep **Fast preview** on by default for the combined view; click
**Plot Full Resolution** when you want ticks, overlays, and the arrays used by
Save Image / Export Script.  See {doc}`../tutorials/gui_web_walkthrough`.

## Palettes and color vision

Full guide (generators, GUI menu map, `plt.show` tips): {doc}`palettes`.

```python
print(mcf.list_palettes())
s.apply_palette('pob')                 # purple / orange / blue, etc.
print(s.colorblind_report())           # flag clashes for common CVD types
# Perceptual auto-N (GUI "Even (auto-N)"):
colors = mcf.suggest_colors(4)
# Hue pattern / classical HSV wheel:
# colors = mcf.colors_for_hue_pattern('triad', n=3)
# colors = mcf.colors_from_hsv(3)

# Scripting preview (pyplot figure — plt.show() works):
res = mcf.preview_palette('pob', n=3, mode='lab', blend='screen',
                          background='black')
# res.show()  or  import matplotlib.pyplot as plt; plt.show()
# res.cvd['ok'], res.colors; save with res.fig.savefig(...)
```

The GUI **Colors** tab applies curated palettes and shows a live color-vision
status.  Side-by-side looks on one target: {doc}`../examples/crab_palette_gallery`.
`preview_palette` is the script form: solid tiles on a chosen background, the
same mode-aware overlapping-circle swatch used on figures, and deuteranopia /
protanopia / tritanopia simulations with ΔE warnings.

## Persistence

```python
s.save_state('session.json')     # paths + settings (not pixels)
s.load_state('session.json')
print(s.export_script())         # standalone recreation script
```

Launch with state:

```python
mcf.gui(state='session.json')
mcf.gui_qt(files=['a.fits', 'b.fits'], colors=['#f00', '#0ff'])
```

Alignment when grids differ: {doc}`wcs_and_alignment`.
Tutorials: {doc}`../tutorials/session_export_script`,
{doc}`../tutorials/save_and_export`.
