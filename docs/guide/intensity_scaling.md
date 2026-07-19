# Intensity scaling

Per-panel stretch and min/max control how each FITS layer is mapped to
greyscale before colorization. This is separate from compositing — see
{doc}`color_compositing` for how tinted layers are merged.

---

## Built-in stretches

`to_grey_rgb` and `rescale_image` accept any name in `mcf.stretch_functions`:

`linear`, `sqrt`, `squared`, `log`, `power`, `sinh`, `asinh`.

The browser and Qt GUIs expose the same list in each panel's stretch dropdown
(signed stretches below are scripting-only for now).  **Zscale** / percentile
sliders set `vmin`/`vmax` the same way as `PanelState.apply_zscale()` and
`min_max=` on `to_grey_rgb`.

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
