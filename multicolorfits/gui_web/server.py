"""
FastAPI backend for the multicolorfits browser GUI.

This is a local, single-user tool: the server holds one McfSession in memory
and the browser page is a thin view onto it.  All endpoints are JSON or PNG.
"""

import io
import os
import tempfile

import base64
import struct

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import image as mpl_image
from matplotlib.backends.backend_agg import FigureCanvasAgg

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles

from ..session import (
    McfSession, STRETCHES, COMBINE_MODES, ALIGN_FRAMES, reproject_available,
    MAX_PANELS, MIN_PANELS,
)
from ..overlays import overlays_available
from ..figures import make_combined_figure
from ..core.colormix import COLORSPACES, BLENDS
from ..core.stack_tools import downsample_for_preview
from ..palettes import PALETTES, list_palettes, list_palette_menu, is_palette_name, rotate_colors_from_base

ALIGN_TARGETS = ('reference',) + tuple(ALIGN_FRAMES)

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
UPLOAD_CACHE_DIR = os.path.join(os.path.expanduser('~'), '.cache', 'multicolorfits', 'uploads')


def _downsample(arr, max_size):
    """Stride-downsample an image so its longest side is <= max_size."""
    if max_size and max(arr.shape[0], arr.shape[1]) > max_size:
        step = int(np.ceil(max(arr.shape[0], arr.shape[1]) / max_size))
        arr = arr[::step, ::step]
    return arr


def _png_response(rgb, max_size=None):
    """Encode a [0..1] RGB array as a PNG response (origin: lower)."""
    rgb = _downsample(np.asarray(rgb), max_size)
    rgb = np.clip(np.nan_to_num(rgb), 0, 1)
    buf = io.BytesIO()
    # Flip vertically: FITS origin is lower-left, PNG origin is upper-left
    mpl_image.imsave(buf, rgb[::-1], format='png')
    return Response(content=buf.getvalue(), media_type='image/png')


def _rgba_png_response(rgba, max_size=None):
    """Encode a [0..1] RGBA array as a PNG response (origin: lower)."""
    rgba = np.asarray(rgba)
    if max_size and max(rgba.shape[:2]) > max_size:
        step = int(np.ceil(max(rgba.shape[:2]) / max_size))
        rgba = rgba[::step, ::step]
    rgba = np.clip(np.nan_to_num(rgba), 0, 1)
    buf = io.BytesIO()
    mpl_image.imsave(buf, rgba[::-1], format='png')
    return Response(content=buf.getvalue(), media_type='image/png')


def create_app(session=None):
    session = session or McfSession()
    app = FastAPI(title='multicolorfits')
    app.state.session = session

    def get_panel(idx: int):
        if idx < 0 or idx >= len(session.panels):
            raise HTTPException(status_code=404, detail='No such panel')
        return session.panels[idx]

    # ---------- state ----------

    @app.get('/api/state')
    def get_state():
        state = session.to_dict()
        state['stretches'] = STRETCHES
        state['overlays_available'] = overlays_available()
        state['max_panels'] = MAX_PANELS
        state['min_panels'] = MIN_PANELS
        return state

    @app.post('/api/panels/add')
    def add_panel(payload: dict = None):
        # Accept an empty body (browser "+ Add panel" posts {}).
        if payload is None:
            payload = {}
        try:
            idx = session.add_panel(color=payload.get('color'),
                                    label=payload.get('label', ''))
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        return {'index': idx, 'state': session.to_dict()}

    @app.post('/api/panels/{idx}/remove')
    def remove_panel(idx: int):
        try:
            session.remove_panel(idx)
        except (ValueError, IndexError) as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        return {'state': session.to_dict()}

    # ---------- grid alignment ----------

    @app.get('/api/grid_status')
    def grid_status():
        report = session.grid_report()
        report['reproject_available'] = reproject_available()
        report['align_targets'] = list(ALIGN_TARGETS)
        return report

    @app.post('/api/align')
    def align(payload: dict):
        target = payload.get('target', 'reference')
        if target not in ALIGN_TARGETS:
            raise HTTPException(status_code=400, detail='target must be one of %s' % (ALIGN_TARGETS,))
        reference = payload.get('reference')
        if not reproject_available():
            # Every align path resamples with reproject; fail early and clearly.
            raise HTTPException(status_code=400,
                                detail="Aligning requires the 'reproject' package. "
                                       'Install with:  pip install multicolorfits[reproject]')
        try:
            result = session.align_panels(target=target, reference=reference)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except ImportError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Alignment failed: %s' % exc)
        return {'result': result, 'state': session.to_dict()}

    # ---------- palettes & colorblind check ----------

    @app.get('/api/palettes')
    def palettes():
        curated = [{'name': name, 'colors': PALETTES[name]} for name in list_palettes()]
        menu = [{'group': g, 'items': [{'label': lbl, 'key': key} for lbl, key in items]}
                for g, items in list_palette_menu()]
        return {'palettes': curated, 'menu': menu}

    @app.post('/api/apply_palette')
    def apply_palette(payload: dict):
        name = payload.get('name', '')
        if not is_palette_name(name):
            raise HTTPException(status_code=400, detail='Unknown palette %r' % name)
        result = session.apply_palette(name)
        report = session.colorblind_report()
        return {'result': result, 'colorblind': report, 'state': session.to_dict()}

    @app.post('/api/rotate_hues')
    def rotate_hues(payload: dict):
        delta = float(payload.get('delta_deg', 0))
        base = payload.get('base_colors')
        result = session.rotate_panel_hues(delta, base_colors=base)
        report = session.colorblind_report()
        return {'result': result, 'colorblind': report, 'state': session.to_dict()}

    @app.post('/api/preview_swatch')
    def preview_swatch(payload: dict):
        """Overlapping-circle swatch preview for hue tuning (does not change panel colors)."""
        delta = float(payload.get('delta_deg', 0))
        base = payload.get('base_colors') or session.panel_colors()
        if not base:
            raise HTTPException(status_code=400, detail='No loaded layers')
        colors = rotate_colors_from_base(list(base), delta)
        idx = session.active_indices()
        labels = [session.panels[i].label or ('Image %d' % (i + 1)) for i in idx]
        size = max(128, min(512, int(getattr(session.compose, 'combo_swatch_size', 320) or 320)))
        data = session.render_combo_swatch_colors(colors, labels=labels, size=size)
        if data is None:
            raise HTTPException(status_code=400, detail='Could not render swatch')
        return _rgba_png_response(data['rgba'])

    @app.get('/api/colorblind')
    def colorblind(min_distance: float = 25.0):
        return session.colorblind_report(min_distance=min_distance)

    # ---------- panel data loading ----------

    @app.post('/api/panel/{idx}/load')
    def load_panel(idx: int, payload: dict):
        panel = get_panel(idx)
        path = os.path.expanduser(payload.get('path', ''))
        if not os.path.isfile(path):
            raise HTTPException(status_code=400, detail='File not found: %s' % path)
        try:
            panel.load_fits(path)
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Could not load FITS: %s' % exc)
        return panel.to_dict()

    @app.post('/api/panel/{idx}/upload')
    async def upload_panel(idx: int, file: UploadFile = File(...)):
        panel = get_panel(idx)
        orig_name = os.path.basename(file.filename or 'upload.fits')
        suffix = os.path.splitext(orig_name)[1] or '.fits'
        if not suffix.lower().startswith('.fit'):
            suffix = '.fits'
        os.makedirs(UPLOAD_CACHE_DIR, exist_ok=True)
        dest = os.path.join(UPLOAD_CACHE_DIR, orig_name)
        try:
            content = await file.read()
            with open(dest, 'wb') as f:
                f.write(content)
            panel.load_fits(dest)
            panel.filepath = dest
        except Exception as exc:
            if os.path.isfile(dest):
                try:
                    os.unlink(dest)
                except OSError:
                    pass
            raise HTTPException(status_code=400, detail='Could not load FITS: %s' % exc)
        return panel.to_dict()

    @app.post('/api/panel/{idx}/clear')
    def clear_panel(idx: int):
        get_panel(idx).clear()
        return {'ok': True}

    # ---------- panel parameters ----------

    @app.post('/api/panel/{idx}/params')
    def set_panel_params(idx: int, payload: dict):
        panel = get_panel(idx)
        if 'stretch' in payload:
            if payload['stretch'] not in STRETCHES:
                raise HTTPException(status_code=400, detail='Unknown stretch')
            panel.stretch = payload['stretch']
        if 'color' in payload:
            panel.color = payload['color']
        if 'label' in payload:
            panel.label = str(payload['label'])
        if 'smooth' in payload:
            panel.smooth = bool(payload['smooth'])
        if 'smooth_sigma' in payload:
            panel.smooth_sigma = float(payload['smooth_sigma'])
        # Absolute limits take precedence over percentiles if both sent
        if 'vmin' in payload or 'vmax' in payload:
            panel.set_limits(payload.get('vmin'), payload.get('vmax'))
        elif 'percent_min' in payload or 'percent_max' in payload:
            panel.set_percentiles(payload.get('percent_min'), payload.get('percent_max'))
        return panel.to_dict()

    @app.post('/api/panel/{idx}/zscale')
    def panel_zscale(idx: int):
        panel = get_panel(idx)
        panel.apply_zscale()
        return panel.to_dict()

    @app.post('/api/panel/{idx}/minmax')
    def panel_minmax(idx: int):
        panel = get_panel(idx)
        panel.reset_minmax()
        return panel.to_dict()

    # ---------- panel views ----------

    @app.get('/api/panel/{idx}/preview.png')
    def panel_preview(idx: int, inverse: bool = False, max_size: int = 0):
        panel = get_panel(idx)
        if not panel.in_use:
            raise HTTPException(status_code=404, detail='No image loaded')
        c = session.compose
        if getattr(c, 'panel_preview_hd', False):
            disp = panel.render_display(gamma=c.gamma, inverse=inverse)
            return _png_response(disp, max_size=max_size or None)
        preview_max = int(getattr(c, 'panel_preview_max_size', 512) or 512)
        disp = panel.render_display(gamma=c.gamma, inverse=inverse,
                                    dtype=np.float32, max_size=preview_max)
        return _png_response(disp, max_size=preview_max)

    @app.get('/api/panel/{idx}/preview_buffer')
    def panel_preview_buffer(idx: int, max_size: int = 512):
        """Downsampled float32 pixel buffer for fast in-browser level previews."""
        panel = get_panel(idx)
        if not panel.in_use or panel.data is None:
            raise HTTPException(status_code=404, detail='No image loaded')
        preview_max = max(64, int(max_size or 512))
        c = session.compose
        if getattr(c, 'panel_preview_hd', False):
            preview_max = max(preview_max, 1024)
        data = downsample_for_preview(panel.data, max_size=preview_max, order=1)
        data = np.asarray(data, dtype=np.float32)
        return {
            'shape': [int(data.shape[0]), int(data.shape[1])],
            'data_min': float(panel.data_min),
            'data_max': float(panel.data_max),
            'bytes': base64.b64encode(data.tobytes()).decode('ascii'),
        }

    @app.get('/api/panel/{idx}/preview_buffer.bin')
    def panel_preview_buffer_bin(idx: int, max_size: int = 768):
        """Binary float32 preview buffer: 8-byte header (ny,nx uint32 LE) + pixels."""
        panel = get_panel(idx)
        if not panel.in_use or panel.data is None:
            raise HTTPException(status_code=404, detail='No image loaded')
        preview_max = max(64, int(max_size or 768))
        c = session.compose
        if getattr(c, 'panel_preview_hd', False):
            preview_max = max(preview_max, 1024)
        data = downsample_for_preview(panel.data, max_size=preview_max, order=1)
        data = np.asarray(data, dtype=np.float32)
        ny, nx = int(data.shape[0]), int(data.shape[1])
        header = struct.pack('<II', ny, nx)
        return Response(content=header + data.tobytes(), media_type='application/octet-stream')

    @app.get('/api/panel/{idx}/histogram')
    def panel_histogram(idx: int, bins: int = 100):
        panel = get_panel(idx)
        if not panel.in_use:
            raise HTTPException(status_code=404, detail='No image loaded')
        finite = panel.data[np.isfinite(panel.data)]
        counts, edges = np.histogram(finite, bins=bins)
        return {'counts': counts.tolist(), 'edges': edges.tolist()}

    @app.get('/api/panel/{idx}/header', response_class=PlainTextResponse)
    def get_header(idx: int):
        return get_panel(idx).header_string()

    @app.post('/api/panel/{idx}/header')
    def set_header(idx: int, payload: dict):
        panel = get_panel(idx)
        try:
            panel.apply_header_string(payload.get('text', ''))
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Invalid header: %s' % exc)
        return {'ok': True}

    # ---------- compose ----------

    @app.post('/api/compose')
    def set_compose(payload: dict):
        from ..session import _json_float
        c = session.compose

        def _f(key, default=None):
            if key not in payload:
                return None
            return _json_float(payload[key], default)

        g = _f('gamma')
        if g is not None: c.gamma = g
        if 'inverse' in payload: c.inverse = bool(payload['inverse'])
        if 'tickcolor' in payload: c.tickcolor = payload['tickcolor']
        if 'facecolor' in payload: c.facecolor = payload['facecolor'] or 'none'
        if 'minorticks' in payload: c.minorticks = bool(payload['minorticks'])
        v = _f('tick_major_size')
        if v is not None: c.tick_major_size = v
        v = _f('tick_minor_size')
        if v is not None: c.tick_minor_size = v
        v = _f('tick_major_width')
        if v is not None: c.tick_major_width = v
        v = _f('tick_minor_width')
        if v is not None: c.tick_minor_width = v
        if 'tick_direction' in payload:
            d = str(payload['tick_direction']).lower()
            if d not in ('in', 'out', 'inout'):
                raise HTTPException(status_code=400, detail='tick_direction must be in, out, or inout')
            c.tick_direction = d
        if 'coord_style' in payload: c.coord_style = payload['coord_style']
        if 'x_format' in payload: c.x_format = payload['x_format']
        if 'y_format' in payload: c.y_format = payload['y_format']
        if 'title' in payload: c.title = payload['title']
        if 'xlabel' in payload: c.xlabel = payload['xlabel']
        if 'ylabel' in payload: c.ylabel = payload['ylabel']
        if 'bare_plot' in payload: c.bare_plot = bool(payload['bare_plot'])
        if 'panel_preview_hd' in payload: c.panel_preview_hd = bool(payload['panel_preview_hd'])
        v = _f('panel_preview_max_size')
        if v is not None: c.panel_preview_max_size = max(64, int(v))
        if 'show_legend' in payload: c.show_legend = bool(payload['show_legend'])
        if 'legend_loc' in payload: c.legend_loc = str(payload['legend_loc'])
        if 'show_combo_swatch' in payload: c.show_combo_swatch = bool(payload['show_combo_swatch'])
        if 'combo_swatch_loc' in payload: c.combo_swatch_loc = str(payload['combo_swatch_loc'])
        if 'combo_swatch_labels' in payload: c.combo_swatch_labels = bool(payload['combo_swatch_labels'])
        v = _f('combo_swatch_label_offset')
        if v is not None: c.combo_swatch_label_offset = v
        v = _f('combo_swatch_inset_scale')
        if v is not None: c.combo_swatch_inset_scale = float(max(0.08, min(0.45, v)))
        v = _f('combo_swatch_size')
        if v is not None: c.combo_swatch_size = max(64, min(1024, int(v)))
        if 'show_band_labels' in payload: c.show_band_labels = bool(payload['show_band_labels'])
        if 'band_labels_loc' in payload: c.band_labels_loc = str(payload['band_labels_loc'])
        if 'show_compass' in payload: c.show_compass = bool(payload['show_compass'])
        if 'compass_loc' in payload: c.compass_loc = str(payload['compass_loc'])
        if 'show_beam' in payload: c.show_beam = bool(payload['show_beam'])
        if 'beam_loc' in payload: c.beam_loc = str(payload['beam_loc'])
        if 'beam_style' in payload: c.beam_style = str(payload['beam_style'])
        if 'show_scale_bar' in payload: c.show_scale_bar = bool(payload['show_scale_bar'])
        v = _f('scale_bar_asec')
        if v is not None: c.scale_bar_asec = v
        v = _f('scale_bar_loc')
        if v is not None: c.scale_bar_loc = int(v)
        if 'scale_bar_color' in payload: c.scale_bar_color = payload['scale_bar_color']
        if 'scale_bar_stroke_color' in payload: c.scale_bar_stroke_color = payload['scale_bar_stroke_color']
        v = _f('scale_bar_stroke_lw')
        if v is not None: c.scale_bar_stroke_lw = v
        if 'combine_mode' in payload:
            if payload['combine_mode'] not in COMBINE_MODES:
                raise HTTPException(status_code=400, detail='combine_mode must be one of %s' % (COMBINE_MODES,))
            c.combine_mode = payload['combine_mode']
        if 'combine_blend' in payload:
            if payload['combine_blend'] not in BLENDS:
                raise HTTPException(status_code=400, detail='combine_blend must be one of %s' % (BLENDS,))
            c.combine_blend = payload['combine_blend']
        if 'combine_background' in payload:
            bg = str(payload['combine_background'] or 'black').strip()
            if bg.lower() not in ('black', 'white', 'transparent', 'none'):
                from matplotlib.colors import to_rgb
                try:
                    to_rgb(bg)
                except (ValueError, TypeError):
                    raise HTTPException(status_code=400,
                                        detail='combine_background must be black/white/transparent or a valid color')
            c.combine_background = bg
        return c.to_dict()

    @app.get('/api/combined.png')
    def combined_png(inverse: bool = False, max_size: int = 1024, preview: bool = False):
        """
        Fast pixel-only combined preview (no WCS axes).

        When preview=True, the compose runs in float32 with the source data
        downsampled to <= max_size pixels first, for quick interactive updates
        on large images.  The full-resolution WCS plot (combined_plot.png) and
        FITS export always use full float64 precision.
        """
        try:
            if preview:
                combined = session.render_combined(inverse=inverse or None,
                                                   preview=True, max_size=max_size)
                # Source already downsampled; no need to shrink the PNG again.
                return _png_response(combined, max_size=None)
            combined = session.render_combined(inverse=inverse or None)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        return _png_response(combined, max_size=max_size)

    @app.get('/api/combined_plot.png')
    def combined_plot_png(dpi: int = 100):
        """Full matplotlib WCS-projected plot of the combined image."""
        try:
            fig = make_combined_figure(session)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except Exception as exc:
            raise HTTPException(status_code=500, detail='Plot failed: %s' % exc)
        buf = io.BytesIO()
        try:
            canvas = FigureCanvasAgg(fig)
            canvas.draw()
            canvas.print_png(buf)
        except Exception as exc:
            raise HTTPException(status_code=500, detail='Plot failed: %s' % exc)
        return Response(content=buf.getvalue(), media_type='image/png')

    # ---------- outputs ----------

    @app.post('/api/save/image')
    def save_image(payload: dict):
        path = os.path.expanduser(payload.get('path', ''))
        if not path:
            raise HTTPException(status_code=400, detail='No save path given')
        dpi = int(payload.get('dpi', 300))
        try:
            fig = make_combined_figure(session)
            transparent = str(session.compose.facecolor).strip().lower() in ('none', 'transparent', '')
            fig.savefig(path, dpi=dpi, bbox_inches='tight',
                        facecolor=fig.get_facecolor(), transparent=transparent)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Could not save: %s' % exc)
        return {'ok': True, 'path': path}

    def _cutout_kwargs_from_payload(payload):
        """Parse shared transparent-cutout options from a JSON body."""
        size = payload.get('size', None)
        if size in ('', None):
            size = None
        else:
            size = int(size)
        soft = float(payload.get('soft_edge', 0) or 0)
        pad = float(payload.get('pad', 0.05) if payload.get('pad') is not None else 0.05)
        crop = payload.get('crop', 'auto')
        if crop in (False, 'false', 'none', 'off'):
            crop = 'none'
        return dict(
            size=size,
            alpha_lo=float(payload.get('alpha_lo', 55)),
            alpha_hi=float(payload.get('alpha_hi', 99.3)),
            alpha_gamma=float(payload.get('alpha_gamma', 0.5)),
            alpha_source=str(payload.get('alpha_source', 'luma') or 'luma'),
            crop=crop,
            pad=pad,
            matte=str(payload.get('matte', 'none') or 'none'),
            soft_edge=soft,
            invert=bool(payload.get('invert', False)),
        )

    @app.post('/api/export/cutout')
    def export_cutout(payload: dict):
        """
        Build a transparent-sky cutout from the current session.

        Body fields: size, alpha_lo/hi/gamma, alpha_source, crop, pad, matte,
        soft_edge, invert, path (optional server path), download (bool).
        If *download* is true (or *path* is empty), returns a PNG body.
        """
        if not session.active_panels():
            raise HTTPException(status_code=400, detail='No fits file loaded yet')
        kwargs = _cutout_kwargs_from_payload(payload or {})
        path = os.path.expanduser((payload or {}).get('path', '') or '')
        download = bool((payload or {}).get('download', not path))
        try:
            rgba = session.export_transparent_cutout(
                savepath=path if path and not download else None, **kwargs)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Cutout failed: %s' % exc)
        if path and not download:
            return {'ok': True, 'path': path,
                    'shape': list(rgba.shape)}
        # PNG download / preview bytes
        buf = io.BytesIO()
        try:
            from PIL import Image
            arr8 = (np.clip(rgba, 0, 1) * 255).astype(np.uint8)
            # Flip to origin=upper for typical image viewers.
            Image.fromarray(arr8[::-1], mode='RGBA').save(buf, format='PNG')
        except ImportError:
            from matplotlib import image as mpl_image
            mpl_image.imsave(buf, np.clip(rgba[::-1], 0, 1), format='png')
        return Response(
            content=buf.getvalue(),
            media_type='image/png',
            headers={'Content-Disposition': 'attachment; filename="multicolorfits_cutout.png"'},
        )

    @app.get('/api/fits_targets')
    def fits_targets():
        return {'targets': session.fits_save_targets()}

    @app.post('/api/save/fits')
    def save_fits(payload: dict):
        path = os.path.expanduser(payload.get('path', ''))
        if not path:
            raise HTTPException(status_code=400, detail='No save path given')
        kind = payload.get('kind', 'combined')
        index = int(payload.get('index', -1))
        overwrite = bool(payload.get('overwrite', True))
        try:
            session.save_fits_target(kind, index=index, savepath=path,
                                     overwrite=overwrite)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc))
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Could not save: %s' % exc)
        return {'ok': True, 'path': path}

    @app.get('/api/params', response_class=PlainTextResponse)
    def get_params():
        return session.params_text()

    @app.get('/api/export_script', response_class=PlainTextResponse)
    def export_script():
        return session.export_script()

    # ---------- session save / restore (JSON) ----------

    @app.get('/api/session.json')
    def session_json():
        import json
        body = json.dumps(session.to_state_dict(), indent=2)
        return Response(
            content=body,
            media_type='application/json',
            headers={'Content-Disposition': 'attachment; filename="multicolorfits_session.json"'},
        )

    @app.post('/api/session/save')
    def session_save(payload: dict):
        path = os.path.expanduser(payload.get('path', ''))
        if not path:
            raise HTTPException(status_code=400, detail='No save path given')
        try:
            session.save_state(path)
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise OSError('Session file was not written')
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Could not save session: %s' % exc)
        return {'ok': True, 'path': path}

    @app.post('/api/session/load')
    def session_load(payload: dict):
        try:
            base_dir = payload.get('base_dir')
            if 'state' in payload:
                warnings = session.load_state(
                    payload['state'], base_dir=base_dir,
                    extra_search_dirs=[UPLOAD_CACHE_DIR])
            else:
                path = os.path.expanduser(payload.get('path', ''))
                if not path:
                    raise HTTPException(status_code=400, detail='No path or state given')
                warnings = session.load_state(
                    path, base_dir=base_dir, extra_search_dirs=[UPLOAD_CACHE_DIR])
        except HTTPException:
            raise
        except Exception as exc:
            raise HTTPException(status_code=400, detail='Could not load session: %s' % exc)
        return {'warnings': warnings, 'state': session.to_dict()}

    @app.post('/api/session/reset')
    def session_reset():
        session.reset()
        return {'state': session.to_dict()}

    @app.get('/api/cursor')
    def cursor_info(x: int, y: int, img_w: int = 0, img_h: int = 0):
        """Sample layer values at a data pixel (optionally map from display coords)."""
        if img_w > 0 and img_h > 0:
            dx, dy = session.display_to_data_pixel(x, y, img_w, img_h)
            if dx is None:
                raise HTTPException(status_code=400, detail='No image loaded')
            x, y = dx, dy
        return session.sample_at_pixel(x, y)

    # ---------- server-side file browsing (local tool) ----------

    @app.get('/api/browse')
    def browse(path: str = '.'):
        # '.' / '' => process CWD (handy default for local GUI use).
        raw = (path or '.').strip() or '.'
        if raw in ('.', './'):
            full = os.getcwd()
        else:
            full = os.path.abspath(os.path.expanduser(raw))
        if not os.path.isdir(full):
            raise HTTPException(status_code=400, detail='Not a directory: %s' % full)
        dirs, files = [], []
        try:
            for entry in sorted(os.listdir(full)):
                if entry.startswith('.'):
                    continue
                p = os.path.join(full, entry)
                if os.path.isdir(p):
                    dirs.append(entry)
                elif entry.lower().endswith(('.fits', '.fit', '.fts', '.fits.gz')):
                    files.append(entry)
        except PermissionError:
            raise HTTPException(status_code=403, detail='Permission denied')
        return {
            'path': full,
            'parent': os.path.dirname(full),
            'dirs': dirs,
            'files': files,
            'cwd': os.getcwd(),
        }

    # ---------- frontend ----------

    @app.middleware('http')
    async def _static_no_cache(request, call_next):
        """Avoid stale panel-preview JS after upgrades (browsers cache /static aggressively)."""
        response = await call_next(request)
        path = request.url.path
        if path == '/' or path.startswith('/static/'):
            if path.endswith(('.js', '.css', '.html')) or path == '/':
                response.headers['Cache-Control'] = 'no-cache, must-revalidate'
        return response

    @app.get('/')
    def index():
        return FileResponse(os.path.join(STATIC_DIR, 'index.html'))

    app.mount('/static', StaticFiles(directory=STATIC_DIR), name='static')

    return app
