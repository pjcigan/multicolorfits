# Intensity scaling

Per-panel stretch and min/max control how each FITS layer is mapped to
greyscale before colorization. This is separate from compositing — see
{doc}`color_compositing` for how tinted layers are merged. Job snippets:
{doc}`../capabilities/index`. Recipes: `mcf.recipes('suggest')`,
`mcf.recipes('levels')`.

---

## Built-in stretches

`to_grey_rgb` and `rescale_image` accept any name in `mcf.stretch_functions`:

`linear`, `sqrt`, `squared`, `log`, `power`, `sinh`, `asinh`.

The browser and Qt GUIs expose the same list in each panel's stretch dropdown
(signed stretches below are scripting-only for now).  **Zscale** sets limits
only. **Auto levels** also picks a stretch — see
{ref}`suggested-levels` below. Percentile sliders set `vmin`/`vmax` the same
way as `min_max=` on `to_grey_rgb`.

---

(suggested-levels)=
## Suggested levels

{func}`multicolorfits.suggest_levels` looks at the pixel distribution and
returns a stretch plus absolute `vmin`/`vmax` to use as the **first** values
on a panel. It is a starting point. It does not lock those fields, and it
does not claim the image is well displayed — especially not an object whose
faint structure and bright peaks span many decades. That case needs more
than one mapping; this function only gets each layer of a composite into a
reasonable first state so you can nudge by eye.

{func}`multicolorfits.describe_image` returns a header summary (when you pass
a header) and the suggestion, and prints both by default. Pass
``verbose=False`` to keep the dict without printing — handy in a loop.
``PanelState.apply_suggested_levels`` writes the stretch and limits onto
one panel. ``McfSession.apply_suggested_levels`` does that for every loaded
panel (the script form of **Auto levels**). The same button in the browser
and Qt GUIs updates one panel and puts the one-line reason in the status
bar. Neither runs on load.

For several layers at once, {func}`multicolorfits.describe_images` runs the
same characterization, suggests colors for the set, and returns parallel
lists ready for ``load_files`` / ``set_display`` (quiet by default).

```python
import multicolorfits as mcf

rec = mcf.suggest_levels(data)
print(rec['stretch'], rec['vmin'], rec['vmax'])
print(rec['reason'])

grey = mcf.to_grey_rgb(
    data, rescalefn=rec['stretch'], scaletype='abs',
    min_max=[rec['vmin'], rec['vmax']])

# Header plus levels (printed):
mcf.describe_image(data, hdr, name='layer')
# Collect without printing:
info = mcf.describe_image(data, hdr, name='layer', verbose=False)
suggested = info['levels']   # stretch, vmin, vmax, reason, …

# Session: Auto levels on every loaded panel (prints the reason if verbose).
recs = s.apply_suggested_levels(verbose=True)

# Values you already chose (scalar broadcasts; sequences match loaded panels):
s.set_display(stretches=['asinh', 'sqrt'], vmins=[0.2, 0.1], vmaxs=[12, 8])
# load_files accepts the same optional stretches / vmins / vmaxs.
```

### Several layers at once

{func}`multicolorfits.describe_images` is the batch form of
{func}`multicolorfits.describe_image`. Pass a mapping
``name -> (data, header)`` (or a sequence of arrays / FITS paths), and it
returns suggested colors plus stretch/limit lists you can feed straight
into a session:

```python
# aligned is name -> (data, header), e.g. after prep_layers / crop
report = mcf.describe_images(aligned)
# report['names'], ['colors'], ['stretches'], ['vmins'], ['vmaxs']
# report['layers'][name] -> {'header', 'levels', 'color'}

# Keep your own colors, or a named palette instead of suggest_colors:
# report = mcf.describe_images(aligned, colors=['#3B5BDB', '#0CA678', '#E03131'])
# report = mcf.describe_images(aligned, palette='pob')

s = mcf.McfSession(n_panels=len(report['names']))
s.load_files(
    paths,   # same order as report['names']
    colors=report['colors'],
    labels=report['names'],
    stretches=report['stretches'],
    vmins=report['vmins'],
    vmaxs=report['vmaxs'],
)
```

Default is quiet (``verbose=False``). Pass ``verbose=True`` to print each
layer's header and levels dump. For a single image, prefer
``describe_image(..., verbose=False)`` over rolling your own loop.

Exact zeros are dropped by default so chip gaps and blanked borders do not
become the floor. Pass ``ignore_zeros=False`` when zero is a real
measurement. An image that is nothing but exact zeros then has nothing left
to measure and raises ``ValueError``. Non-finite pixels are always ignored.
At most 250000 pixels are sampled (fixed seed).

The stretch is only ever ``linear``, ``sqrt``, or ``asinh``. ``log`` is never
chosen: it collapses at or below zero. If the image crosses zero, the reason
mentions ``symlog`` as a script option, but it is not selected (the GUI
stretch list does not include it). The optional ``symmetric_log`` stretch is
never selected.

Rules, applied in order:

1. Fewer than 16 usable pixels: ``linear`` between the sample min and max.
2. At least 1% of the sample is negative, or the background sits on the
   floor: ``asinh``. Floor means the 5th percentile is at or below zero, or
   the lower half is piled up (the median is less than 3× the 5th
   percentile) while the 99.5th percentile is at least 20× the median. A
   distribution that fills many decades is not a floor — the lower half is
   spread out — and uses the next rule instead.
3. Otherwise the span is ``log10(p99.5 / p5)``. Under about 1 decade:
   ``linear``. Under about 2: ``sqrt``. Wider: ``asinh``.
4. Limits are a mild clip: ``vmin`` is the 1st percentile, ``vmax`` is the
   99.9th. If the 99.9th is at least twice the 99.5th, the bright tail is
   thin, ``vmax`` stays at the 99.5th, and the reason says the core will
   saturate.
5. If the span is more than about 4 decades, the stretch is unchanged and
   the reason notes that one stretch will not show faint structure and
   bright peaks together.

A constant image still returns ``linear``, with a tiny range so a later
stretch does not divide by zero.

**Zscale** remains the mid-range clip (IRAF-style). It does not choose a
stretch, and it is the wrong first guess when a thin bright tail or a floor
full of faint structure is what you care about. Auto levels does not replace
it; both buttons stay on the panel.

---

## Signed flux through zero

Background-subtracted maps often span negative noise through positive signal.
`log` cannot handle negatives; `asinh` works but treats the near-zero region
differently than symmetric logarithmic scaling.

| `rescalefn` | Dependency | Behaviour |
|-------------|------------|-----------|
| `'asinh'` | core | Smooth; already in GUI |
| `'symlog'` | core (matplotlib `SymLogNorm`) | Piecewise linear + log through zero |
| `'symmetric_log'` | optional `[pysymlog]` | C¹-continuous symlog (smoother near zero) |

```python
import multicolorfits as mcf

grey = mcf.to_grey_rgb(
    data,
    rescalefn='symlog',
    scaletype='abs',
    min_max=[vmin, vmax],
)

# For histograms / residual side panels
norm = mcf.make_norm('symmetric_log', vmin=-0.5, vmax=2.0)  # needs pysymlog
```

Install the optional extra:

```bash
pip install multicolorfits[pysymlog]
```

Colormap-heavy display (diverging maps, grey wedges for shallow negatives) is
better handled in **skyplothelper** or a dedicated colormap workflow — mcf
colorizes single-hue layers, not scalar colormap images.  For a minimal
matplotlib diverging map over a signed stretch, use ``TwoSlopeNorm`` (or
``CenteredNorm``) with a diverging colormap on the greyscale data before (or
instead of) tinting.
