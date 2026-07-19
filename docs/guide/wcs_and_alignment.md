# WCS, headers, and alignment

Colorize-then-combine only works when every layer shares the **same pixel
grid**. Headers and WCS decide that grid; crop / reproject / align tools get
you there. API detail: {doc}`../api/wcs_fits`.

---

## `McfHeader`

Each `PanelState` stores its FITS header as an
{class}`~multicolorfits.McfHeader` — a thin wrapper with cached WCS, celestial
frame, and pixel scale, plus `.matches_grid(other)` for cheap same-shape /
same-WCS checks.

```python
import multicolorfits as mcf

h = mcf.McfHeader(header)          # or mcf.as_mcfheader(header)
h.wcs                             # astropy WCS (cached)
h.matches_grid(other_header)
h.minimal_wcs_header()            # portable WCS-only sidecar
mcf.save_header(h, 'layer.hdr')   # / load_header(...)
```

`load_fits` / `PanelState.load_fits` squeeze stray Stokes / singleton axes so
2D imaging layers land as `(ny, nx)` without a manual `squeeze_image` in most
cases.

---

## Same grid vs needs alignment

| Situation | What to do |
|-----------|------------|
| Already registered (same shape + WCS) | Load and combine |
| Different pointing / pixel scale / rotation | Reproject onto one reference |
| Want a common celestial frame (e.g. Galactic north up) | Align with a frame target |

In either GUI a mismatch warning offers **Align layers…** (needs
`[reproject]`). Scripts use the session helper or the stack API:

```python
# Session — preserves stretch / color / labels on each panel
s.align_panels(target='reference')          # onto first loaded panel
# s.align_panels(target='galactic')         # common grid, Galactic north up

# Arrays + headers already in hand
aligned = mcf.align_stack(
    [(data_a, hdr_a), (data_b, hdr_b), (data_c, hdr_c)],
    reference=0,
)
# aligned = mcf.align_stack(..., frame='icrs')  # optional frame
```

Single-layer helpers: `reproject_image`, `reproject_cube`,
`reproject_to_frame`, `reproject_to_galactic`, `optimal_common_header`.

Worked narrative: {doc}`../examples/wlm` and tutorial
{doc}`../tutorials/wlm_align_reproject`.

---

## Crop before (or instead of) a full reproject

Trim to a shared sky or pixel window when companions are larger than the
region you care about:

```python
data_c, hdr_c = mcf.crop_image(data, header, xbounds=[100, 400], ybounds=[120, 420])
data_c, hdr_c = mcf.crop_image_sky(data, header, ra_bounds=..., dec_bounds=...)
```

Cubes: `crop_cube` / `crop_cube_sky`. After an interactive GUI crop, **Save
Session** and reload, or crop in script then `load_files` / `set_data`.

---

## Coordinates and overlays

Pixel ↔ sky: `sky_to_pixel`, `pixel_to_sky`, plus sexagesimal helpers
(`sex2dec`, …). Combined figures can draw a compass, beam, and scale bar when
the `[overlays]` extra (skyplothelper) is installed — see
{doc}`figures_and_mosaics` and {doc}`../api/overlays`.
