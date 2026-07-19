"""Tests for combined-figure construction, including canvas (facecolor) control."""

import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.figures import (
    contrast_color, make_combined_figure, add_band_labels, add_swatch_inset,
)

from conftest import make_test_header, make_test_data


@pytest.fixture
def loaded_session():
    s = mcf.McfSession()
    s.panels[0].set_data(make_test_data(), make_test_header())
    s.panels[0].color = '#FF0000'
    return s


@pytest.fixture
def multi_session():
    s = mcf.McfSession()
    for i, col in enumerate(['#FF0000', '#00FF00', '#0000FF']):
        s.panels[i].set_data(make_test_data(seed=i), make_test_header())
        s.panels[i].color = col
        s.panels[i].label = 'band%d' % i
    return s


class TestAnnotations:
    def test_figure_with_swatch_and_band_labels(self, multi_session):
        multi_session.compose.show_combo_swatch = True
        multi_session.compose.show_band_labels = True
        multi_session.compose.combo_swatch_labels = True
        fig = make_combined_figure(multi_session)
        assert fig is not None

    def test_add_band_labels_creates_texts(self, loaded_session):
        fig = make_combined_figure(loaded_session)
        ax = fig.axes[0]
        before = len(ax.texts)
        add_band_labels(ax, [('#FF0000', 'Ha'), ('#00FF00', 'OIII')], loc='upper left')
        assert len(ax.texts) >= before + 2

    def test_add_swatch_inset_adds_axes(self, multi_session):
        fig = make_combined_figure(multi_session)
        ax = fig.axes[0]
        n_before = len(ax.child_axes)
        add_swatch_inset(ax, multi_session.render_combo_swatch(size=64))
        assert len(ax.child_axes) == n_before + 1

    def test_add_swatch_inset_none_noop(self, loaded_session):
        fig = make_combined_figure(loaded_session)
        ax = fig.axes[0]
        n_before = len(ax.child_axes)
        add_swatch_inset(ax, None)
        assert len(ax.child_axes) == n_before


class TestContrastColor:
    def test_light_backgrounds_get_black_text(self):
        assert contrast_color('white') == 'black'
        assert contrast_color('#ffffff') == 'black'
        assert contrast_color('0.9') == 'black'

    def test_dark_backgrounds_get_white_text(self):
        assert contrast_color('black') == 'white'
        assert contrast_color('#001f3f') == 'white'
        assert contrast_color('0.1') == 'white'

    def test_transparent_defaults_to_black(self):
        assert contrast_color('none') == 'black'
        assert contrast_color('') == 'black'

    def test_unparseable_defaults_to_black(self):
        assert contrast_color('not-a-color') == 'black'


class TestFacecolorRendering:
    def test_default_facecolor_is_white(self):
        assert mcf.ComposeState().facecolor == 'white'

    def test_white_canvas(self, loaded_session):
        loaded_session.compose.facecolor = 'white'
        fig = make_combined_figure(loaded_session)
        assert fig.get_facecolor() == pytest.approx((1.0, 1.0, 1.0, 1.0))

    def test_colored_canvas(self, loaded_session):
        loaded_session.compose.facecolor = '#204060'
        fig = make_combined_figure(loaded_session)
        r, g, b, a = fig.get_facecolor()
        assert (round(r, 2), round(g, 2), round(b, 2), a) == (0.13, 0.25, 0.38, 1.0)

    def test_transparent_canvas_has_zero_alpha(self, loaded_session):
        loaded_session.compose.facecolor = 'none'
        fig = make_combined_figure(loaded_session)
        assert fig.get_facecolor()[3] == 0.0

    def test_transparent_canvas_axes_patch_transparent(self, loaded_session):
        loaded_session.compose.facecolor = 'none'
        fig = make_combined_figure(loaded_session)
        assert fig.axes[0].patch.get_alpha() == 0.0


class TestBarePlot:
    def test_bare_plot_hides_frame(self, loaded_session):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        loaded_session.compose.bare_plot = True
        fig = make_combined_figure(loaded_session)
        FigureCanvasAgg(fig).draw()
        ax = fig.axes[0]
        assert ax.get_xlabel() == ''
        assert ax.get_ylabel() == ''
        assert ax.get_title() == ''
        assert len(ax.coords[0].ticklabels.get_text()) == 0
        assert ax.get_legend() is None

    def test_apply_bare_plot_style(self):
        c = mcf.ComposeState(show_legend=True, show_compass=True)
        mcf.apply_bare_plot_style(c, transparent=True)
        assert c.bare_plot is True
        assert c.facecolor == 'none'
        assert c.show_legend is False
        assert c.show_compass is False

    def test_session_apply_bare_plot_style(self):
        s = mcf.McfSession()
        s.apply_bare_plot_style(transparent=True)
        assert s.compose.bare_plot is True
        assert s.compose.facecolor == 'none'

    def test_bare_plot_session_roundtrip(self, loaded_session, tmp_path):
        loaded_session.compose.bare_plot = True
        p = tmp_path / 'state.json'
        loaded_session.save_state(str(p))
        s2 = mcf.McfSession()
        s2.load_state(str(p))
        assert s2.compose.bare_plot is True

    def test_export_script_bare_plot(self, loaded_session):
        loaded_session.compose.bare_plot = True
        script = loaded_session.export_script()
        assert 'bare_plot=True' in script


class TestAxisLabelsAndTicks:
    def test_custom_axis_labels(self, loaded_session):
        loaded_session.compose.xlabel = 'Right Ascension'
        loaded_session.compose.ylabel = 'Declination'
        fig = make_combined_figure(loaded_session)
        ax = fig.axes[0]
        assert ax.get_xlabel() == 'Right Ascension'
        assert ax.get_ylabel() == 'Declination'

    def test_auto_axis_labels_from_ctype(self, loaded_session):
        fig = make_combined_figure(loaded_session)
        ax = fig.axes[0]
        assert ax.get_xlabel() == 'RA'
        assert ax.get_ylabel() == 'DEC'

    def test_tick_direction_out(self, loaded_session):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        loaded_session.compose.tick_direction = 'out'
        loaded_session.compose.tick_major_size = 10.0
        fig = make_combined_figure(loaded_session)
        FigureCanvasAgg(fig).draw()
        rap = fig.axes[0].coords[0]
        assert len(rap.ticklabels.get_text()) >= 1

    def test_rotated_wcs_shows_ticklabels(self, loaded_session):
        """Rotated fields (e.g. NGC 602) place ticks on non-default frame edges."""
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        hdr = loaded_session.common_header.astropy.copy()
        hdr['CROTA2'] = 90.0
        loaded_session.panels[0].header = loaded_session.panels[0].header.__class__(hdr)
        loaded_session.compose.xlabel = 'RA'
        loaded_session.compose.ylabel = 'DEC'
        fig = make_combined_figure(loaded_session)
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        ax = fig.axes[0]
        rap, decp = ax.coords[0], ax.coords[1]
        ra_text = [t for side in rap.ticklabels.text for t in rap.ticklabels.text[side]]
        dec_text = [t for side in decp.ticklabels.text for t in decp.ticklabels.text[side]]
        assert any(str(t).strip() for t in ra_text), 'expected RA tick labels on rotated field'
        assert any(str(t).strip() for t in dec_text), 'expected Dec tick labels on rotated field'
        assert ax.coords[0].get_axislabel() == 'RA'
        assert ax.coords[1].get_axislabel() == 'DEC'

    def test_tick_direction_inout_does_not_fail(self, loaded_session):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        loaded_session.compose.tick_direction = 'inout'
        fig = make_combined_figure(loaded_session)
        FigureCanvasAgg(fig).draw()

    def test_compose_tick_fields_roundtrip(self, loaded_session, tmp_path):
        c = loaded_session.compose
        c.xlabel = 'X'
        c.ylabel = 'Y'
        c.tick_major_size = 9.0
        c.tick_minor_size = 3.5
        c.tick_direction = 'out'
        p = tmp_path / 'state.json'
        loaded_session.save_state(str(p))
        s2 = mcf.McfSession()
        s2.load_state(str(p))
        assert s2.compose.xlabel == 'X'
        assert s2.compose.tick_direction == 'out'
        assert s2.compose.tick_major_size == 9.0


class TestFacecolorSerialization:
    def test_to_dict_includes_facecolor(self):
        assert 'facecolor' in mcf.ComposeState().to_dict()

    def test_save_restore_roundtrips_facecolor(self, loaded_session, tmp_path):
        loaded_session.compose.facecolor = '#123456'
        p = tmp_path / 'state.json'
        loaded_session.save_state(str(p))
        s2 = mcf.McfSession()
        s2.load_state(str(p))
        assert s2.compose.facecolor == '#123456'

    def test_export_script_includes_facecolor_and_labelcolor(self, loaded_session):
        loaded_session.compose.facecolor = '#204060'
        line = [l for l in loaded_session.export_script().splitlines()
                if 'plot_combined_rgb' in l][0]
        assert "facecolor='#204060'" in line
        assert "labelcolor='white'" in line

    def test_export_script_transparent(self, loaded_session):
        loaded_session.compose.facecolor = 'none'
        line = [l for l in loaded_session.export_script().splitlines()
                if 'plot_combined_rgb' in l][0]
        assert "facecolor='none'" in line
        assert "labelcolor='black'" in line


@pytest.fixture
def two_layer_session():
    s = mcf.McfSession()
    s.panels[0].set_data(make_test_data(), make_test_header(), filepath='/d/halpha.fits')
    s.panels[0].color = '#FF0000'
    s.panels[1].set_data(make_test_data(), make_test_header(), filepath='/d/oiii.fits')
    s.panels[1].color = '#0000FF'
    s.panels[1].label = '[OIII]'
    return s


class TestChannelLegend:
    def test_default_no_legend(self, two_layer_session):
        fig = make_combined_figure(two_layer_session)
        assert fig.axes[0].get_legend() is None

    def test_legend_enabled(self, two_layer_session):
        two_layer_session.compose.show_legend = True
        fig = make_combined_figure(two_layer_session)
        assert fig.axes[0].get_legend() is not None

    def test_legend_labels(self, two_layer_session):
        two_layer_session.compose.show_legend = True
        fig = make_combined_figure(two_layer_session)
        texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert texts == ['halpha', '[OIII]']  # filename stem + explicit label

    def test_legend_entries_fallback_label(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(), make_test_header())  # no filepath
        s.panels[0].label = ''
        assert s.legend_entries() == [('#FFFFFF', 'Image 1')]

    def test_legend_loc_used(self, two_layer_session):
        two_layer_session.compose.show_legend = True
        two_layer_session.compose.legend_loc = 'lower left'
        fig = make_combined_figure(two_layer_session)
        assert fig.axes[0].get_legend() is not None

    def test_label_default_from_filename(self):
        s = mcf.McfSession()
        s.panels[0].set_data(make_test_data(), make_test_header(), filepath='/x/ngc602_ir.fits')
        assert s.panels[0].label == 'ngc602_ir'

    def test_compose_legend_roundtrip(self, two_layer_session, tmp_path):
        two_layer_session.compose.show_legend = True
        two_layer_session.compose.legend_loc = 'lower right'
        p = tmp_path / 'state.json'
        two_layer_session.save_state(str(p))
        s2 = mcf.McfSession()
        s2.load_state(str(p))
        assert s2.compose.show_legend is True
        assert s2.compose.legend_loc == 'lower right'

    def test_label_roundtrip_with_real_file(self, tmp_path):
        import astropy.io.fits as pyfits
        fp = tmp_path / 'chan.fits'
        pyfits.writeto(str(fp), make_test_data(), make_test_header())
        s = mcf.McfSession()
        s.panels[0].load_fits(str(fp))
        s.panels[0].label = 'Custom Name'
        state = tmp_path / 'state.json'
        s.save_state(str(state))
        s2 = mcf.McfSession()
        s2.load_state(str(state))
        assert s2.panels[0].label == 'Custom Name'

    def test_export_script_legend_live(self, two_layer_session):
        """When show_legend is on, the exported plot call passes a live legend."""
        two_layer_session.compose.show_legend = True
        script = two_layer_session.export_script()
        # legend is emitted as a real kwarg on the plot call (not a comment)
        plot_line = [l for l in script.splitlines()
                     if l.startswith('mcf.plot_combined_rgb')][0]
        assert 'legend=[' in plot_line
        assert 'legend_loc=' in plot_line
        assert 'halpha' in script and '[OIII]' in script
        compile(script, '<x>', 'exec')

    def test_export_script_no_legend_when_off(self, two_layer_session):
        two_layer_session.compose.show_legend = False
        script = two_layer_session.export_script()
        assert 'legend=[' not in script

    def test_export_script_includes_overlays_and_label_offset(self, two_layer_session):
        two_layer_session.compose.show_compass = True
        two_layer_session.compose.show_scale_bar = True
        two_layer_session.compose.show_combo_swatch = True
        two_layer_session.compose.combo_swatch_labels = True
        two_layer_session.compose.combo_swatch_label_offset = 0.55
        script = two_layer_session.export_script()
        assert 'show_compass=True' in script
        assert 'show_scale_bar=True' in script
        assert 'swatch_label_offset=0.55' in script
        assert 'multicolorfits[overlays]' in script

    def test_plot_combined_rgb_legend_param(self, two_layer_session, tmp_path):
        """The scripting plot helper draws a legend when passed one."""
        combined = two_layer_session.render_combined()
        hdr = two_layer_session.common_header.astropy
        out = tmp_path / 'legended.png'
        mcf.plot_combined_rgb(
            combined, hdr, 'legend test', str(out),
            legend=[('#FF0000', 'Halpha'), ('#0000FF', '[OIII]')],
            legend_loc='upper left', dpi=50)
        assert out.exists()


class TestSaveImageTransparency:
    def test_transparent_png_has_alpha_channel(self, loaded_session, tmp_path):
        loaded_session.compose.facecolor = 'none'
        fig = make_combined_figure(loaded_session)
        out = tmp_path / 'transparent.png'
        fig.savefig(str(out), dpi=50, bbox_inches='tight',
                    facecolor=fig.get_facecolor(), transparent=True)
        from matplotlib import image as mpl_image
        arr = mpl_image.imread(str(out))
        assert arr.shape[-1] == 4  # RGBA
        # At least some fully-transparent margin pixels
        assert (arr[..., 3] == 0).any()
