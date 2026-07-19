"""Qt desktop GUI tests -- skipped automatically when PySide6 is unavailable.

Focus: the compose-controls resync so a pre-configured session
(``mcf.gui_qt(session=s)`` / ``MainWindow(session=s)``) is reflected in the
widgets without the sync clobbering the restored state.
"""
import os

import numpy as np
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
pytest.importorskip('PySide6')

import multicolorfits as mcf  # noqa: E402
from multicolorfits.qt_gui.app import MainWindow  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope='module')
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


def _configured_session():
    s = mcf.McfSession()
    d = np.random.default_rng(0).random((16, 16))
    for i, c in enumerate(['#FF3030', '#FFD000', '#2E6BFF']):
        s.panels[i].set_data(d * (i + 1))
        s.panels[i].vmin = 0
        s.panels[i].vmax = 3
        s.panels[i].color = c
    return s


class TestQtComposeResync:
    def test_gui_qt_stays_callable_after_submodule_import(self):
        import multicolorfits as mcf
        from multicolorfits.qt_gui.app import MainWindow  # noqa: F401
        assert callable(mcf.gui_qt)

    def test_controls_reflect_passed_session(self, qapp):
        s = _configured_session()
        s.compose.combine_mode = 'ryb'
        s.compose.combine_background = 'transparent'
        s.compose.gamma = 1.7
        s.compose.inverse = True
        s.compose.title = 'hello'
        s.compose.tickcolor = '#123456'
        w = MainWindow(session=s)
        assert w.combine_mode_combo.currentData() == 'ryb'
        assert w.combine_bg_combo.currentData() == 'transparent'
        assert abs(w.gamma_spin.value() - 1.7) < 1e-6
        assert w.inverse_check.isChecked() is True
        assert w.title_edit.text() == 'hello'
        assert w.tickcolor_edit.text() == '#123456'

    def test_resync_does_not_mutate_state(self, qapp):
        # Syncing RYB + black must NOT auto-flip the background to white
        # (that convenience only applies to a live user selection).
        s = mcf.McfSession()
        s.compose.combine_mode = 'ryb'
        s.compose.combine_background = 'black'
        w = MainWindow(session=s)
        assert s.compose.combine_background == 'black'
        assert w.combine_bg_combo.currentData() == 'black'

    def test_custom_background_shows_swatch(self, qapp):
        s = mcf.McfSession()
        s.compose.combine_background = '#123456'
        w = MainWindow(session=s)
        assert w.combine_bg_combo.currentData() == 'custom'
        # isHidden() reflects the explicit setVisible() state without needing the
        # (never-shown) top-level window to be visible.
        assert w.combine_bg_btn.isHidden() is False

    def test_noncustom_background_hides_swatch(self, qapp):
        s = mcf.McfSession()
        s.compose.combine_background = 'white'
        w = MainWindow(session=s)
        assert w.combine_bg_combo.currentData() == 'white'
        assert w.combine_bg_btn.isHidden() is True

    def test_blend_enabled_only_for_perceptual(self, qapp):
        s = mcf.McfSession(); s.compose.combine_mode = 'lab'
        assert MainWindow(session=s).combine_blend_combo.isEnabled() is True
        s2 = mcf.McfSession(); s2.compose.combine_mode = 'ryb'
        assert MainWindow(session=s2).combine_blend_combo.isEnabled() is False
        s3 = mcf.McfSession()  # default rgb
        assert MainWindow(session=s3).combine_blend_combo.isEnabled() is False

    def test_loaded_panels_reflected(self, qapp):
        s = _configured_session()
        w = MainWindow(session=s)
        assert w.panel_widgets[0].color_edit.text() == '#FF3030'
        assert w.panel_widgets[1].color_edit.text() == '#FFD000'

    def test_live_selection_still_updates(self, qapp):
        # After construction, user interaction should still work normally.
        s = mcf.McfSession()
        w = MainWindow(session=s)
        idx = w.combine_mode_combo.findData('lab')
        w.combine_mode_combo.setCurrentIndex(idx)
        assert s.compose.combine_mode == 'lab'
        assert w.combine_blend_combo.isEnabled() is True

    def test_swatch_and_band_label_controls_reflect_session(self, qapp):
        s = _configured_session()
        s.compose.show_combo_swatch = True
        s.compose.combo_swatch_loc = 'upper left'
        s.compose.combo_swatch_labels = True
        s.compose.show_band_labels = True
        s.compose.band_labels_loc = 'lower right'
        w = MainWindow(session=s)
        assert w.swatch_check.isChecked() is True
        assert w.swatch_loc_combo.currentText() == 'upper left'
        assert w.swatch_labels_check.isChecked() is True
        assert w.band_labels_check.isChecked() is True
        assert w.band_labels_loc_combo.currentText() == 'lower right'

    def test_swatch_toggle_updates_state(self, qapp):
        s = mcf.McfSession()
        w = MainWindow(session=s)
        w.swatch_check.setChecked(True)
        assert s.compose.show_combo_swatch is True
        w.band_labels_check.setChecked(True)
        assert s.compose.show_band_labels is True

    def test_overlay_controls_reflect_session(self, qapp):
        s = _configured_session()
        s.compose.show_compass = True
        s.compose.show_beam = True
        s.compose.show_scale_bar = True
        w = MainWindow(session=s)
        assert w.compass_check.isChecked() is True
        assert w.beam_check.isChecked() is True
        assert w.scale_bar_check.isChecked() is True

    def test_overlay_toggle_updates_state(self, qapp):
        s = mcf.McfSession()
        w = MainWindow(session=s)
        w.compass_check.setChecked(True)
        w.beam_check.setChecked(True)
        w.scale_bar_check.setChecked(True)
        assert s.compose.show_compass is True
        assert s.compose.show_beam is True
        assert s.compose.show_scale_bar is True

    def test_overlay_placement_controls_reflect_session(self, qapp):
        s = _configured_session()
        s.compose.compass_loc = 'upper right'
        s.compose.beam_style = 'ellipse'
        s.compose.scale_bar_asec = 12.0
        s.compose.combo_swatch_label_offset = 0.4
        w = MainWindow(session=s)
        assert w.compass_loc_combo.currentText() == 'upper right'
        assert w.beam_style_combo.currentText() == 'ellipse'
        assert w.scale_bar_asec_spin.value() == pytest.approx(12.0)
        assert w.swatch_offset_spin.value() == pytest.approx(0.4)

    def test_light_theme_uses_light_palette(self, qapp):
        from multicolorfits.qt_gui.app import _light_palette
        from PySide6.QtGui import QPalette
        w = MainWindow(session=mcf.McfSession())
        w._apply_theme('light')
        assert w._theme == 'light'
        pal = qapp.palette()
        assert pal.color(QPalette.Window) == _light_palette().color(QPalette.Window)

    def test_scale_bar_color_controls_reflect_session(self, qapp):
        s = mcf.McfSession()
        s.compose.scale_bar_color = '#ffcc00'
        s.compose.scale_bar_stroke_color = '#ffffff'
        s.compose.scale_bar_stroke_lw = 2.5
        w = MainWindow(session=s)
        assert w.scale_bar_color_edit.text() == '#ffcc00'
        assert w.scale_bar_stroke_edit.text() == '#ffffff'
        assert w.scale_bar_stroke_lw_spin.value() == pytest.approx(2.5)

    def test_add_panel_tab_grows_widgets(self, qapp):
        w = MainWindow(session=mcf.McfSession())
        assert w.tabs.count() == 4
        w.add_panel_tab()
        assert w.tabs.count() == 5
        assert len(w.session.panels) == 5
        assert w.tabs.tabText(4) == 'Image 5'
        w.remove_panel_tab(4)
        assert w.tabs.count() == 4
        assert len(w.session.panels) == 4

    def test_mainwindow_honors_extra_session_panels(self, qapp):
        s = mcf.McfSession(n_panels=6)
        w = MainWindow(session=s)
        assert w.tabs.count() == 6
        assert len(w.panel_widgets) == 6
