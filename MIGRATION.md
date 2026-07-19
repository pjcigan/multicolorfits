# Migrating from multicolorfits v2.1.x to v3

## TL;DR

* **Scripting code runs unchanged.**  Every v2.x public name still resolves
  with the same signature and numerical behavior (`greyRGBize_image`,
  `colorize_image`, `combine_multicolor`, `hexinv`, `saveRGBfits`,
  `makesimpleheader`, `reproject2D`, `cropfits2D`, ...).
* **Many functions have clearer new names** (e.g. `to_grey_rgb`,
  `hex_complement`, `save_rgb_fits`, `reproject_image`).  The old names are
  kept as silent, fully functional aliases — see **Renamed functions** below.
  New code should use the new names; the aliases will emit `DeprecationWarning`
  in a later 3.x release and be removed no earlier than v4.0.
* **The GUI stack is completely replaced.**  traits, traitsui, pyface, and
  PyQt are no longer dependencies (and are never imported).  The primary GUI
  is now browser-based; a PySide6 desktop GUI is available as a second option.
* **GUI dependencies are opt-in extras** rather than hard requirements:
  `pip install multicolorfits` gives you a headless scripting library.
  Most users will want ``pip install "multicolorfits[all]"`` (both GUIs,
  alignment, overlays, and related helpers).

## What changed

### Dependencies

| v2.1.x (always required)          | v3                                   |
|-----------------------------------|---------------------------------------|
| numpy, astropy, scipy, matplotlib, scikit-image | same (only these)      |
| traits, traitsui, pyface          | removed                               |
| PyQt4 / PyQt5 / PyQt6             | removed                               |
| --                                | `[web]` extra: fastapi, uvicorn       |
| --                                | `[qt]` extra: PySide6                 |
| `[reproject]` extra               | same (now a lazy import; no import-time warning) |

This resolves the long-standing v2.x issues with TraitsUI on Python >= 3.11
and the PyQt6 segfaults, since none of those packages are used anymore.

**Conda note.**  v3 no longer depends on PyQt, but the optional `[qt]` extra
installs **PySide6** (a pip Qt6 stack).  Pip-installing Qt6 into a conda env
that already has conda-forge Qt5 / PyQt5 can break interactive
`matplotlib` `plt.show()` (heap corruption once the GUI event loop starts),
even when Agg/headless plots still work.  Prefer `[web]` in shared conda
Qt5 environments; put `[qt]` in a dedicated env.  Full audit / repair steps
are in the Sphinx installation guide (*Conda environments and the `[qt]` /
PySide6 extra*), also available as
[`docs/installation.md`](https://github.com/pjcigan/multicolorfits/blob/master/docs/installation.md)
in the repository.

### Launching the GUI

```python
# v2.1.x
import multicolorfits as mcf
mcf.mcf_gui()

# v3 -- easiest: pip install "multicolorfits[all]", then:
import multicolorfits as mcf
mcf.gui()          # browser GUI (also needs [web] if not using [all])
# mcf.gui_qt()     # desktop GUI (also needs [qt] if not using [all])
# mcf.mcf_gui()    # still works as an alias for gui()
```

Console scripts are also installed: `mcf-web` and `mcf-qt`.

### Module structure

v2.x was a single 2,400-line `multicolorfits.py`.  v3 is a package with the
computation split into `multicolorfits.core.*` submodules, but everything is
re-exported at the top level, so `mcf.<anything>` from v2.x still resolves.

### Renamed functions (v2 name -> v3 name)

Most public functions were renamed in v3.0 for clarity and PEP 8 consistency.
Every old name remains available as a fully functional alias (no warnings
emitted yet); the full old -> new mapping is importable as
`multicolorfits.LEGACY_ALIASES`, and all aliases live in one module
(`multicolorfits/compat.py`) for easy future removal.

| v2.x name | v3 name |
|---|---|
| `greyRGBize_image` | `to_grey_rgb` |
| `hexinv` | `hex_complement` |
| `rgb_to_hsv_vectorized` / `hsv_to_rgb_vectorized` | `rgb_to_hsv` / `hsv_to_rgb` |
| `rgb_to_hsl_vectorized` / `hsl_to_rgb_vectorized` | `rgb_to_hsl` / `hsl_to_rgb` |
| `rgb_to_ryb_vectorized` / `ryb_to_rgb_vectorized` | `rgb_to_ryb` / `ryb_to_rgb` |
| `rgb_to_cmyk_vectorized` / `cmyk_to_rgb_vectorized` | `rgb_to_cmyk` / `cmyk_to_rgb` |
| `image_value` | `layer_coverage` |
| `plotsinglemulticolorRGB` | `plot_combined_rgb` |
| `comparemulticolorRGB_pureRGB` | `compare_multicolor_vs_rgb` |
| `saveRGBfits` | `save_rgb_fits` |
| `scaling_fns` | `stretch_functions` |
| `drawProgressBar` | `draw_progress_bar` |
| `nanpercofscore` | `nan_percentile_of_score` |
| `makesimpleheader` | `make_simple_header` |
| `getcdelts` / `getcdmatrix` | `get_cdelts` / `get_cd_matrix` |
| `getdegperpix` / `getasecperpix` / `getsteradperpix` | `deg_per_pixel` / `arcsec_per_pixel` / `sterad_per_pixel` |
| `convsky2pix` / `convpix2sky` | `sky_to_pixel` / `pixel_to_sky` |
| `angulardistance` | `angular_distance` |
| `beampars_asec_fromhdr` | `beam_params_arcsec` |
| `pixperbeam_from_hdr` | `pixels_per_beam` |
| `force_hdr_to_2D` / `force_hdr_to_3D` | `force_header_2d` / `force_header_3d` |
| `force_hdr_floats` | `force_header_floats` |
| `cropfits2D` / `cropfits3D` | `crop_image` / `crop_cube` |
| `cropfits2D_coords` / `cropfits3D_coords` | `crop_image_sky` / `crop_cube_sky` |
| `reproject2D` / `reproject3D` | `reproject_image` / `reproject_cube` |

Names that were already clear (`colorize_image`, `combine_multicolor`,
`hex_to_rgb`, `deg2dms`, `smooth_image`, ...) are unchanged.

### Behavior notes

* Importing `multicolorfits` no longer prints warnings when `reproject` or
  `kapteyn` are missing.  Instead, `reproject2D()`/`reproject3D()` raise a
  clear `ImportError` with install instructions only if you actually call
  them without the needed package.
* `drawProgressBar` no longer raises `NameError` (a missing `sys` import in
  v2.x).
* `combine_multicolor` no longer emits divide-by-zero warnings (and garbage
  output) in the corner case where all channels saturate in inverse mode.
* `smooth_image` works with modern scikit-image (`channel_axis` replaced the
  removed `multichannel` keyword; both are supported).
* `dec2sex` no longer depends on `np.NZERO` (removed in numpy 2.0).

### New (additive) API

These did not exist in v2.x and are safe to ignore:

* `mcf.rescale_image(data, stretch, vmin, vmax)` -- the core scaling step as
  a standalone function.
* `mcf.zscale_limits(data)` -- zscale (vmin, vmax) without the GUI.
* `mcf.PanelState`, `mcf.ComposeState`, `mcf.McfSession` -- the state model /
  controller that both GUIs share; usable as a mid-level scripting API,
  including `export_script()`, `save_state()` / `load_state()`, and `load_files()`.

## Known gaps vs v2.1.x GUI

* The per-panel matplotlib zoom toolbar of the old GUI is replaced by
  simple previews in the panels (browser GUI); the combined plot in the Qt
  GUI keeps the full matplotlib navigation toolbar.
* **Cursor readout** on the combined view shows per-layer flux at `(x, y)` plus
  RA/Dec when WCS is available (web status line + Qt status bar).

## GUI session files and scripted launch

v3 sessions are plain **JSON** files (`.json`), portable between machines as long
as the FITS paths still resolve:

```python
import multicolorfits as mcf

s = mcf.McfSession()
s.load_files(['img1.fits', 'img2.fits'], colors=['#C11B17', '#4CC417'])
s.compose.combine_mode = 'lab'
s.save_state('my_session.json')

# Later — restore in GUI or script
mcf.gui(state='my_session.json')
mcf.gui_qt(state='my_session.json')

# Or preload FITS paths directly (no JSON file)
mcf.gui(files=['a.fits', 'b.fits'], colors=['#f00', '#0f0'])
```

The `mcf-web` CLI accepts the same: `mcf-web --state my_session.json` or
`mcf-web --files a.fits b.fits --colors '#f00' '#0f0'`.
