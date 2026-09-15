"""
PySide6 desktop GUI for multicolorfits.

A deliberately thin wrapper around the shared McfSession controller
(multicolorfits/session.py): all image processing lives there, this module
only builds widgets and forwards events.  No traits, no pyface.
"""

import os
import sys

import numpy as np

from PySide6.QtCore import Qt, QTimer, QSettings
from PySide6.QtGui import QColor, QImage, QPainter, QPalette, QPixmap
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QColorDialog, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox, QFileDialog, QFormLayout, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QSlider,
    QSpinBox, QSplitter, QStatusBar, QTabWidget, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from ..session import (
    McfSession, STRETCHES, ALIGN_FRAMES, reproject_available,
    MAX_PANELS, MIN_PANELS,
)
from ..figures import setup_combined_axes, refresh_wcs_ticklabels, axis_labels_for_compose
from ..palettes import PALETTES, list_palette_menu, colors_for_hue_pattern
from .hue_wheel import HueWheelDialog

PREVIEW_MAX = 512  # downsample panel previews to at most this many pixels/side


def _dark_palette():
    """Warm dark QPalette aligned with the browser GUI dark theme."""
    bg = QColor('#2c2a28')
    base = QColor('#232120')
    fg = QColor('#e8e4de')
    dim = QColor('#a39e96')
    accent = QColor('#4c9be8')
    p = QPalette()
    p.setColor(QPalette.Window, bg)
    p.setColor(QPalette.WindowText, fg)
    p.setColor(QPalette.Base, base)
    p.setColor(QPalette.AlternateBase, QColor('#353230'))
    p.setColor(QPalette.ToolTipBase, base)
    p.setColor(QPalette.ToolTipText, fg)
    p.setColor(QPalette.Text, fg)
    p.setColor(QPalette.Button, QColor('#353230'))
    p.setColor(QPalette.ButtonText, fg)
    p.setColor(QPalette.BrightText, QColor('#ff6b6b'))
    p.setColor(QPalette.Highlight, accent)
    p.setColor(QPalette.HighlightedText, QColor('#10131a'))
    p.setColor(QPalette.Link, accent)
    p.setColor(QPalette.PlaceholderText, dim)
    for role in (QPalette.Text, QPalette.ButtonText, QPalette.WindowText):
        p.setColor(QPalette.Disabled, role, dim)
    return p


def _light_palette():
    """Light QPalette aligned with the browser GUI light theme."""
    bg = QColor('#f1f3f7')
    base = QColor('#ffffff')
    fg = QColor('#1b202a')
    dim = QColor('#5c6674')
    accent = QColor('#2f7fd0')
    p = QPalette()
    p.setColor(QPalette.Window, bg)
    p.setColor(QPalette.WindowText, fg)
    p.setColor(QPalette.Base, base)
    p.setColor(QPalette.AlternateBase, QColor('#e6eaf0'))
    p.setColor(QPalette.ToolTipBase, base)
    p.setColor(QPalette.ToolTipText, fg)
    p.setColor(QPalette.Text, fg)
    p.setColor(QPalette.Button, QColor('#e6eaf0'))
    p.setColor(QPalette.ButtonText, fg)
    p.setColor(QPalette.BrightText, QColor('#c0392b'))
    p.setColor(QPalette.Highlight, accent)
    p.setColor(QPalette.HighlightedText, QColor('#ffffff'))
    p.setColor(QPalette.Link, accent)
    p.setColor(QPalette.PlaceholderText, dim)
    for role in (QPalette.Text, QPalette.ButtonText, QPalette.WindowText):
        p.setColor(QPalette.Disabled, role, dim)
    return p


def _array_to_qpixmap(rgb):
    """[0..1] RGB float array (origin lower) -> QPixmap."""
    arr = np.clip(np.nan_to_num(rgb), 0, 1)
    if max(arr.shape[:2]) > PREVIEW_MAX:
        step = int(np.ceil(max(arr.shape[:2]) / PREVIEW_MAX))
        arr = arr[::step, ::step]
    arr8 = (arr[::-1] * 255).astype(np.uint8)  # flip: FITS origin is lower-left
    h, w, _ = arr8.shape
    arr8 = np.ascontiguousarray(arr8)
    img = QImage(arr8.data, w, h, 3 * w, QImage.Format_RGB888)
    return QPixmap.fromImage(img.copy())


class HeaderDialog(QDialog):
    """Plain-text FITS header editor."""

    def __init__(self, parent, header_text):
        super().__init__(parent)
        self.setWindowTitle('FITS Header')
        self.resize(700, 550)
        layout = QVBoxLayout(self)
        self.editor = QPlainTextEdit()
        self.editor.setPlainText(header_text)
        self.editor.setStyleSheet('font-family: monospace;')
        layout.addWidget(self.editor)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def text(self):
        return self.editor.toPlainText()


class TextViewDialog(QDialog):
    """Read-only text display (params / exported script)."""

    def __init__(self, parent, title, text, save_suffix=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(700, 550)
        layout = QVBoxLayout(self)
        editor = QPlainTextEdit()
        editor.setPlainText(text)
        editor.setReadOnly(True)
        editor.setStyleSheet('font-family: monospace;')
        layout.addWidget(editor)
        row = QHBoxLayout()
        if save_suffix:
            save_btn = QPushButton('Save to file…')

            def do_save():
                path, _ = QFileDialog.getSaveFileName(self, 'Save', '', '*%s' % save_suffix)
                if path:
                    with open(path, 'w') as f:
                        f.write(text)
            save_btn.clicked.connect(do_save)
            row.addWidget(save_btn)
        row.addStretch()
        close_btn = QPushButton('Close')
        close_btn.clicked.connect(self.accept)
        row.addWidget(close_btn)
        layout.addLayout(row)


class PathDialog(QDialog):
    """Load/save path dialog with pasteable line edit + Browse… (like panel FITS)."""

    def __init__(self, parent, title, *, start_dir='', filters='JSON (*.json);;All files (*)',
                 save=False, default_name=''):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._save = save
        self._filters = filters
        self._start_dir = start_dir or ''
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel('Path (paste or type, or Browse…)'))
        row = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText('/path/to/session.json')
        if default_name and start_dir:
            self.path_edit.setText(os.path.join(start_dir, default_name))
        elif default_name:
            self.path_edit.setText(default_name)
        row.addWidget(self.path_edit)
        browse_btn = QPushButton('Browse…')
        browse_btn.clicked.connect(self._browse)
        row.addWidget(browse_btn)
        layout.addLayout(row)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.path_edit.setFocus()
        self.path_edit.selectAll()

    def _browse(self):
        start = self.path_edit.text().strip() or self._start_dir
        if start and not os.path.isdir(start):
            start = os.path.dirname(start) or self._start_dir
        if self._save:
            path, _ = QFileDialog.getSaveFileName(self, self.windowTitle(), start, self._filters)
        else:
            path, _ = QFileDialog.getOpenFileName(self, self.windowTitle(), start, self._filters)
        if path:
            self.path_edit.setText(path)

    def path(self):
        return self.path_edit.text().strip()


class PanelWidget(QWidget):
    """Controls + preview for one image layer; edits one PanelState."""

    def __init__(self, main, index):
        super().__init__()
        self.main = main
        self.index = index
        self._updating = False  # guard against feedback loops while syncing UI
        self._perc_timer = QTimer(self)
        self._perc_timer.setSingleShot(True)
        self._perc_timer.setInterval(200)
        self._perc_timer.timeout.connect(self.on_percentiles)

        layout = QVBoxLayout(self)

        # --- file row ---
        file_row = QHBoxLayout()
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText('/path/to/image.fits')
        self.file_edit.returnPressed.connect(self.load_from_edit)
        browse_btn = QPushButton('Browse…')
        browse_btn.clicked.connect(self.browse)
        file_row.addWidget(self.file_edit)
        file_row.addWidget(browse_btn)
        layout.addLayout(file_row)

        # --- preview ---
        self.preview = QLabel('no image')
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setMinimumSize(300, 300)
        self.preview.setStyleSheet('background: black; color: #888;')
        self.preview.setScaledContents(False)
        layout.addWidget(self.preview, stretch=1)

        # --- scale & color ---
        grid = QGridLayout()
        grid.addWidget(QLabel('Scale'), 0, 0)
        self.stretch_combo = QComboBox()
        self.stretch_combo.addItems(STRETCHES)
        self.stretch_combo.currentTextChanged.connect(self.on_stretch)
        grid.addWidget(self.stretch_combo, 0, 1)

        grid.addWidget(QLabel('Color'), 0, 2)
        color_row = QHBoxLayout()
        self.color_btn = QPushButton()
        self.color_btn.setFixedWidth(36)
        self.color_btn.clicked.connect(self.pick_color)
        self.color_edit = QLineEdit('#FFFFFF')
        self.color_edit.setFixedWidth(80)
        self.color_edit.editingFinished.connect(self.on_color_text)
        color_row.addWidget(self.color_btn)
        color_row.addWidget(self.color_edit)
        grid.addLayout(color_row, 0, 3)

        grid.addWidget(QLabel('Min'), 1, 0)
        self.vmin_edit = QLineEdit()
        self.vmin_edit.editingFinished.connect(self.on_limits)
        grid.addWidget(self.vmin_edit, 1, 1)
        grid.addWidget(QLabel('Max'), 1, 2)
        self.vmax_edit = QLineEdit()
        self.vmax_edit.editingFinished.connect(self.on_limits)
        grid.addWidget(self.vmax_edit, 1, 3)

        grid.addWidget(QLabel('Label'), 2, 0)
        self.label_edit = QLineEdit()
        self.label_edit.setPlaceholderText('(filename)')
        self.label_edit.setToolTip('Channel name shown in the combined-figure legend')
        self.label_edit.editingFinished.connect(self.on_label)
        grid.addWidget(self.label_edit, 2, 1, 1, 3)
        layout.addLayout(grid)

        # --- percentile sliders (x100 for slider integer resolution) ---
        for name, label in [('pmin', 'Min %'), ('pmax', 'Max %')]:
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            slider = QSlider(Qt.Horizontal)
            slider.setRange(0, 10000)
            slider.setValue(0 if name == 'pmin' else 10000)
            val_label = QLabel('0.00' if name == 'pmin' else '100.00')
            val_label.setFixedWidth(52)
            slider.valueChanged.connect(lambda v, lab=val_label: lab.setText('%.2f' % (v / 100.)))
            slider.valueChanged.connect(self._schedule_percentiles)
            slider.sliderReleased.connect(self._flush_percentiles)
            row.addWidget(slider)
            row.addWidget(val_label)
            layout.addLayout(row)
            setattr(self, name + '_slider', slider)

        # --- buttons ---
        btn_row = QHBoxLayout()
        for label, slot, tip in [
            ('Min/Max', self.on_minmax, 'Reset limits to full data min/max'),
            ('Zscale', self.on_zscale, 'Set limits with the zscale algorithm'),
            ('Auto levels', self.on_auto_levels,
             'Suggest a stretch and vmin/vmax from the pixel distribution. '
             'Starting point only; does not lock the fields.'),
            ('Header', self.edit_header, 'View / edit the FITS header'),
            ('Clear', self.on_clear, 'Clear this panel'),
            ('Remove', self.on_remove, 'Remove this panel tab'),
        ]:
            b = QPushButton(label)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            btn_row.addWidget(b)
            if label == 'Remove':
                self.remove_btn = b
        layout.addLayout(btn_row)

        # --- smoothing ---
        smooth_row = QHBoxLayout()
        self.smooth_check = QCheckBox('Smooth')
        self.smooth_check.toggled.connect(self.on_smooth)
        self.smooth_sigma = QDoubleSpinBox()
        self.smooth_sigma.setRange(0.1, 100.)
        self.smooth_sigma.setValue(3.0)
        self.smooth_sigma.setSingleStep(0.5)
        self.smooth_sigma.valueChanged.connect(self.on_smooth)
        smooth_row.addWidget(self.smooth_check)
        smooth_row.addWidget(QLabel('\u03c3 (pixels)'))
        smooth_row.addWidget(self.smooth_sigma)
        smooth_row.addStretch()
        layout.addLayout(smooth_row)

        self._set_color_swatch('#FFFFFF')

    # ---------- helpers ----------

    @property
    def panel(self):
        return self.main.session.panels[self.index]

    def _set_color_swatch(self, hexcolor):
        self.color_btn.setStyleSheet('background-color: %s;' % hexcolor)

    def sync_from_state(self):
        """Reflect the PanelState in the widgets (without re-triggering slots)."""
        p = self.panel
        self._updating = True
        try:
            self.file_edit.setText(p.filepath or '')
            if not p.in_use:
                self.preview.setPixmap(QPixmap())
                self.preview.setText('no image')
            self.stretch_combo.setCurrentText(p.stretch)
            self.vmin_edit.setText('%.6g' % p.vmin)
            self.vmax_edit.setText('%.6g' % p.vmax)
            self.pmin_slider.setValue(int(round(p.percent_min * 100)))
            self.pmax_slider.setValue(int(round(p.percent_max * 100)))
            self.color_edit.setText(p.color)
            self._set_color_swatch(p.color)
            self.label_edit.setText(p.label)
            self.smooth_check.setChecked(p.smooth)
            self.smooth_sigma.setValue(p.smooth_sigma)
        finally:
            self._updating = False

    def refresh_preview(self):
        p = self.panel
        if not p.in_use:
            self.preview.setPixmap(QPixmap())
            self.preview.setText('no image')
            return
        try:
            c = self.main.session.compose
            if getattr(c, 'panel_preview_hd', False):
                disp = p.render_display(gamma=c.gamma)
            else:
                preview_max = int(getattr(c, 'panel_preview_max_size', PREVIEW_MAX) or PREVIEW_MAX)
                disp = p.render_display(gamma=c.gamma,
                                        dtype=np.float32, max_size=preview_max)
        except Exception as exc:
            self.main.show_status('Preview failed: %s' % exc)
            return
        pix = _array_to_qpixmap(disp)
        self.preview.setText('')
        self.preview.setPixmap(pix.scaled(self.preview.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.panel.in_use:
            self.refresh_preview()

    # ---------- slots ----------

    def browse(self):
        start = ''
        if self.panel.filepath:
            start = os.path.dirname(self.panel.filepath)
        elif self.main is not None:
            start = self.main.last_file_dir()
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open FITS file', start,
            'FITS files (*.fits *.fit *.fts *.fits.gz);;All files (*)')
        if path:
            self.file_edit.setText(path)
            if self.main is not None:
                self.main.remember_file_dir(path)
            self.load(path)

    def load_from_edit(self):
        path = self.file_edit.text().strip()
        if path:
            self.load(path)

    def load(self, path):
        try:
            self.panel.load_fits(path)
        except Exception as exc:
            QMessageBox.warning(self, 'Load failed', 'Could not load FITS file:\n%s' % exc)
            return
        self.sync_from_state()
        self.refresh_preview()
        self.main.refresh_grid_status()
        self.main.refresh_cvd()
        self.main.show_status('Image %d: loaded %s  %s' % (self.index + 1, path, self.panel.data.shape))

    def on_stretch(self, stretch):
        if self._updating: return
        self.panel.stretch = stretch
        self.maybe_refresh()

    def pick_color(self):
        color = QColorDialog.getColor(QColor(self.panel.color), self, 'Image color')
        if color.isValid():
            self.panel.color = color.name().upper()
            self.color_edit.setText(self.panel.color)
            self._set_color_swatch(self.panel.color)
            self.main.refresh_cvd()
            self.maybe_refresh()

    def on_color_text(self):
        if self._updating: return
        v = self.color_edit.text().strip()
        if not v.startswith('#'):
            v = '#' + v
        if len(v) == 7:
            self.panel.color = v.upper()
            self._set_color_swatch(v)
            self.main.refresh_cvd()
            self.maybe_refresh()
        else:
            self.main.show_status('Color must be a #RRGGBB hex string')

    def on_label(self):
        if self._updating: return
        self.panel.label = self.label_edit.text().strip()
        c = self.main.session.compose
        if (c.show_legend or c.show_band_labels or c.show_combo_swatch
                or c.show_compass or c.show_beam or c.show_scale_bar):
            self.main._replot_axes_if_shown()

    def on_limits(self):
        if self._updating or not self.panel.in_use: return
        try:
            vmin = float(self.vmin_edit.text())
            vmax = float(self.vmax_edit.text())
        except ValueError:
            return
        self.panel.set_limits(vmin, vmax)
        self.sync_from_state()
        self.maybe_refresh()

    def on_percentiles(self):
        if self._updating or not self.panel.in_use: return
        self.panel.set_percentiles(self.pmin_slider.value() / 100., self.pmax_slider.value() / 100.)
        self.sync_from_state()
        self.maybe_refresh()

    def _schedule_percentiles(self):
        if self._updating or not self.panel.in_use: return
        self._perc_timer.start()

    def _flush_percentiles(self):
        self._perc_timer.stop()
        self.on_percentiles()

    def on_minmax(self):
        if not self.panel.in_use: return
        self.panel.reset_minmax()
        self.sync_from_state()
        self.maybe_refresh()
        self.main.show_status('Scale reset to min/max')

    def on_zscale(self):
        if not self.panel.in_use: return
        self.panel.apply_zscale()
        self.sync_from_state()
        self.maybe_refresh()
        self.main.show_status('Min/max determined by zscale')

    def on_auto_levels(self):
        if not self.panel.in_use:
            return
        try:
            rec = self.panel.apply_suggested_levels()
        except ValueError as exc:
            QMessageBox.warning(self, 'Auto levels', str(exc))
            return
        self.sync_from_state()
        self.maybe_refresh()
        self.main.show_status(rec.get('reason') or 'Suggested stretch and limits')

    def on_smooth(self, *_):
        if self._updating: return
        self.panel.smooth = self.smooth_check.isChecked()
        self.panel.smooth_sigma = self.smooth_sigma.value()
        if self.panel.in_use:
            self.maybe_refresh()

    def edit_header(self):
        dlg = HeaderDialog(self, self.panel.header_string())
        hdr = self.panel.header
        if hdr is not None:
            bits = []
            if hdr.shape:
                bits.append('%d\u00d7%d' % (hdr.shape[1], hdr.shape[0]))
            if hdr.is_celestial:
                bits.append(hdr.frame)
                if hdr.pixscale_asec:
                    bits.append('%.3g"/px' % hdr.pixscale_asec)
            if bits:
                dlg.setWindowTitle('FITS Header  (%s)' % ', '.join(bits))
        if dlg.exec() == QDialog.Accepted:
            try:
                self.panel.apply_header_string(dlg.text())
                self.main.show_status('Header updated')
            except Exception as exc:
                QMessageBox.warning(self, 'Invalid header', str(exc))

    def on_clear(self):
        self.panel.clear()
        self.file_edit.clear()
        self.sync_from_state()
        self.refresh_preview()
        self.main.refresh_grid_status()
        self.main.refresh_cvd()
        self.main.show_status('Image %d cleared' % (self.index + 1))

    def on_remove(self):
        self.main.remove_panel_tab(self.index)

    def maybe_refresh(self):
        if not self.panel.in_use:
            return
        if self.main.autorefresh_check.isChecked():
            self.refresh_preview()
        # Keep the combined live preview in sync (no-op unless it's enabled)
        self.main.schedule_live_preview()


class MainWindow(QMainWindow):

    def __init__(self, session=None):
        super().__init__()
        self.session = session if session is not None else McfSession()
        self.setWindowTitle('MultiColorFits')
        self.resize(1350, 850)

        central = QWidget()
        self.setCentralWidget(central)
        outer = QVBoxLayout(central)

        # --- top bar ---
        top = QHBoxLayout()
        self.autorefresh_check = QCheckBox('Auto-refresh previews')
        self.autorefresh_check.setChecked(True)
        top.addWidget(self.autorefresh_check)
        self.panel_preview_hd_check = QCheckBox('Full-res panel previews')
        self.panel_preview_hd_check.setToolTip(
            'Panel thumbnails use full FITS resolution (slower when adjusting sliders). '
            'Default is a fast downsampled float32 preview.')
        self.panel_preview_hd_check.toggled.connect(self.on_panel_preview_hd)
        top.addWidget(self.panel_preview_hd_check)
        self.live_preview_check = QCheckBox('Fast preview')
        self.live_preview_check.setToolTip(
            'Fast, downsampled float32 preview (panels + combined) that updates as '
            'you adjust settings. Click "Plot Full Resolution" or uncheck for the '
            'full-resolution WCS plot.')
        self.live_preview_check.toggled.connect(self.on_live_preview_toggled)
        top.addWidget(self.live_preview_check)

        top.addStretch()
        self.theme_btn = QPushButton()
        self.theme_btn.setFixedWidth(40)
        self.theme_btn.clicked.connect(self.toggle_theme)
        top.addWidget(self.theme_btn)
        outer.addLayout(top)

        # --- grid-mismatch warning bar (hidden unless layers need aligning) ---
        self.grid_bar = QFrame()
        self.grid_bar.setObjectName('gridBar')
        self.grid_bar.setStyleSheet(
            '#gridBar { background: #6b4a12; border: 1px solid #e8a04c; border-radius: 5px; }'
            '#gridBar QLabel { color: #ffe9c7; }')
        gb = QHBoxLayout(self.grid_bar)
        gb.setContentsMargins(10, 5, 10, 5)
        self.grid_label = QLabel()
        self.grid_label.setWordWrap(True)
        gb.addWidget(QLabel('\u26a0'))
        gb.addWidget(self.grid_label, stretch=1)
        self.align_btn = QPushButton('Align layers\u2026')
        self.align_btn.clicked.connect(self.open_align_dialog)
        gb.addWidget(self.align_btn)
        self.grid_bar.setVisible(False)
        outer.addWidget(self.grid_bar)

        splitter = QSplitter(Qt.Horizontal)
        outer.addWidget(splitter, stretch=1)

        # --- left: panel tabs + add/remove ---
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        self.tabs = QTabWidget()
        self.panel_widgets = []
        self._build_panel_tabs()
        left_layout.addWidget(self.tabs, stretch=1)
        panel_btns = QHBoxLayout()
        self.add_panel_btn = QPushButton('+ Add panel')
        self.add_panel_btn.setToolTip('Add another image panel (up to %d)' % MAX_PANELS)
        self.add_panel_btn.clicked.connect(self.add_panel_tab)
        panel_btns.addWidget(self.add_panel_btn)
        panel_btns.addStretch()
        left_layout.addLayout(panel_btns)
        splitter.addWidget(left)

        # --- right: combined figure + controls ---
        right = QWidget()
        rlayout = QVBoxLayout(right)

        self.compose_tabs = self._build_compose_tabs()
        rlayout.addWidget(self.compose_tabs)

        self.figure = Figure(figsize=(7, 7))
        self.canvas = FigureCanvas(self.figure)
        self.mpl_toolbar = NavigationToolbar2QT(self.canvas, right)
        rlayout.addWidget(self.canvas, stretch=1)
        rlayout.addWidget(self.mpl_toolbar)

        btns = QHBoxLayout()
        plot_btn = QPushButton('Plot Full Resolution')
        plot_btn.setDefault(True)
        plot_btn.clicked.connect(self.plot_combined)
        btns.addWidget(plot_btn)
        for label, slot in [
            ('Load Session…', self.load_session),
            ('Save Session…', self.save_session),
            ('Save Image…', self.save_image),
            ('Save FITS…', self.save_fits),
            ('Export transparent…', self.export_transparent_cutout),
            ('Reset…', self.reset_session),
            ('Params', self.show_params),
            ('Export Script', self.export_script),
        ]:
            b = QPushButton(label)
            b.clicked.connect(slot)
            btns.addWidget(b)
        rlayout.addLayout(btns)

        splitter.addWidget(right)
        splitter.setSizes([420, 900])

        self.setStatusBar(QStatusBar())
        self._cursor_readout = QLabel('')
        self.statusBar().addPermanentWidget(self._cursor_readout)
        self._motion_cid = None
        self._combined_plot_ready = False
        self._draw_placeholder()

        # Debounce timer for the fast live preview
        self._live_timer = QTimer(self)
        self._live_timer.setSingleShot(True)
        self._live_timer.timeout.connect(self._live_preview)

        # Apply the saved (or default) light/dark theme
        self._theme = None
        saved = QSettings('multicolorfits', 'gui_qt').value('theme', 'dark')
        self._apply_theme(saved if saved in ('dark', 'light') else 'dark')

        # Reflect the (possibly pre-configured) session in all controls.
        self.sync_from_session()
        self.refresh_grid_status()
        self._update_palette_swatch()
        self.refresh_cvd()

    LIVE_PREVIEW_MAX = 1024  # downsample source to this many pixels/side for live preview

    # ---------- recent file directories ----------

    def last_file_dir(self):
        """Directory for the next file dialog (last used, else CWD)."""
        settings = QSettings('multicolorfits', 'gui_qt')
        d = settings.value('last_file_dir', '')
        if d and os.path.isdir(d):
            return d
        return os.getcwd()

    def remember_file_dir(self, path):
        if not path:
            return
        d = path if os.path.isdir(path) else os.path.dirname(path)
        if d and os.path.isdir(d):
            QSettings('multicolorfits', 'gui_qt').setValue('last_file_dir', d)

    # ---------- theme ----------

    def _apply_theme(self, theme):
        app = QApplication.instance()
        if app is not None:
            if app.style().objectName().lower() != 'fusion':
                app.setStyle('Fusion')
            app.setPalette(_dark_palette() if theme == 'dark' else _light_palette())
        self._theme = theme
        canvas_bg = '#232120' if theme == 'dark' else '#f1f3f7'
        self.canvas.setStyleSheet('background: %s;' % canvas_bg)
        self.mpl_toolbar.setStyleSheet('')
        # Show the icon of the theme you'd switch TO.
        self.theme_btn.setText('\u2600' if theme == 'dark' else '\u263e')
        self.theme_btn.setToolTip('Switch to %s theme'
                                  % ('light' if theme == 'dark' else 'dark'))
        try:
            QSettings('multicolorfits', 'gui_qt').setValue('theme', theme)
        except Exception:
            pass

    def toggle_theme(self):
        self._apply_theme('light' if self._theme == 'dark' else 'dark')

    # ---------- panel tabs ----------

    def _build_panel_tabs(self, select_index=None):
        """Rebuild left-hand Image N tabs to match ``session.panels``."""
        current = self.tabs.currentIndex() if self.tabs.count() else 0
        while self.tabs.count():
            w = self.tabs.widget(0)
            self.tabs.removeTab(0)
            if w is not None:
                w.deleteLater()
        self.panel_widgets = []
        for i in range(len(self.session.panels)):
            pw = PanelWidget(self, i)
            self.panel_widgets.append(pw)
            self.tabs.addTab(pw, 'Image %d' % (i + 1))
        if select_index is None:
            select_index = current
        if self.tabs.count():
            self.tabs.setCurrentIndex(max(0, min(int(select_index), self.tabs.count() - 1)))
        self._update_panel_tab_actions()

    def _update_panel_tab_actions(self):
        n = len(self.session.panels)
        if hasattr(self, 'add_panel_btn'):
            self.add_panel_btn.setEnabled(n < MAX_PANELS)
        can_remove = n > MIN_PANELS
        for pw in self.panel_widgets:
            if hasattr(pw, 'remove_btn'):
                pw.remove_btn.setEnabled(can_remove)

    def add_panel_tab(self):
        try:
            idx = self.session.add_panel()
        except ValueError as exc:
            QMessageBox.warning(self, 'Add panel', str(exc))
            return
        self._build_panel_tabs(select_index=idx)
        self.show_status('Added image panel %d' % (idx + 1))

    def remove_panel_tab(self, index):
        if len(self.session.panels) <= MIN_PANELS:
            QMessageBox.information(self, 'Remove panel',
                                    'Need at least one image panel.')
            return
        panel = self.session.panels[index]
        if panel.in_use:
            reply = QMessageBox.question(
                self, 'Remove panel',
                'Remove Image %d and unload its FITS data?\n'
                'Later panels will renumber.' % (index + 1),
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return
        try:
            self.session.remove_panel(index)
        except (ValueError, IndexError) as exc:
            QMessageBox.warning(self, 'Remove panel', str(exc))
            return
        self._build_panel_tabs(select_index=max(0, index - 1))
        self.refresh_grid_status()
        self.refresh_cvd()
        self.schedule_live_preview()
        self.show_status('Removed panel %d' % (index + 1))

    # ---------- state sync ----------

    def sync_from_session(self):
        """Reflect the full session state (panels + compose) in every widget.

        Used at startup so a pre-configured session (e.g. ``mcf.gui_qt(session=s)``)
        is shown correctly, and reusable after any future restore action.
        """
        if len(self.panel_widgets) != len(self.session.panels):
            self._build_panel_tabs()
        self._sync_compose_controls()
        for pw in self.panel_widgets:
            pw.sync_from_state()
            if pw.panel.in_use:
                pw.refresh_preview()
        self._update_panel_tab_actions()

    def _sync_compose_controls(self):
        """Reflect session.compose in the compose widgets without re-triggering
        their slots (which would otherwise clobber restored values)."""
        c = self.session.compose
        widgets = [self.gamma_spin, self.inverse_check, self.combine_mode_combo,
                   self.combine_blend_combo, self.combine_bg_combo, self.coord_combo,
                   self.minorticks_check, self.tickcolor_edit, self.facecolor_edit,
                   self.title_edit, self.xlabel_edit, self.ylabel_edit,
                   self.bare_plot_check, self.panel_preview_hd_check,
                   self.tick_major_size_spin, self.tick_minor_size_spin,
                   self.tick_major_width_spin, self.tick_minor_width_spin,
                   self.tick_direction_combo,
                   self.legend_check, self.legend_loc_combo,
                   self.swatch_check, self.swatch_loc_combo, self.swatch_labels_check,
                   self.band_labels_check, self.band_labels_loc_combo,
                   self.swatch_offset_spin,
                   self.swatch_inset_scale_spin, self.swatch_size_spin,
                   self.compass_check, self.beam_check, self.scale_bar_check,
                   self.compass_loc_combo, self.beam_loc_combo, self.beam_style_combo,
                   self.scale_bar_asec_spin, self.scale_bar_loc_combo,
                   self.scale_bar_color_edit, self.scale_bar_stroke_edit,
                   self.scale_bar_stroke_lw_spin]
        for wdg in widgets:
            wdg.blockSignals(True)
        try:
            self.gamma_spin.setValue(float(c.gamma))
            self.inverse_check.setChecked(bool(c.inverse))
            i = self.combine_mode_combo.findData(c.combine_mode)
            if i >= 0:
                self.combine_mode_combo.setCurrentIndex(i)
            self.combine_blend_combo.setCurrentText(c.combine_blend)
            bg_key = str(c.combine_background or 'black').strip().lower()
            if bg_key in ('black', 'white', 'transparent', 'none'):
                j = self.combine_bg_combo.findData('transparent' if bg_key == 'none' else bg_key)
            else:
                j = self.combine_bg_combo.findData('custom')
            if j >= 0:
                self.combine_bg_combo.setCurrentIndex(j)
            self.coord_combo.setCurrentText(
                'Sexagesimal' if str(c.coord_style).startswith('sex') else 'Decimal')
            self.minorticks_check.setChecked(bool(c.minorticks))
            self.tickcolor_edit.setText(str(c.tickcolor))
            self.facecolor_edit.setText(str(c.facecolor))
            self.title_edit.setText(c.title or '')
            self.bare_plot_check.setChecked(bool(getattr(c, 'bare_plot', False)))
            self.panel_preview_hd_check.setChecked(bool(getattr(c, 'panel_preview_hd', False)))
            self.xlabel_edit.setText(c.xlabel or '')
            self.ylabel_edit.setText(c.ylabel or '')
            self.tick_major_size_spin.setValue(float(getattr(c, 'tick_major_size', 8.0)))
            self.tick_minor_size_spin.setValue(float(getattr(c, 'tick_minor_size', 4.0)))
            self.tick_major_width_spin.setValue(float(getattr(c, 'tick_major_width', 1.0)))
            self.tick_minor_width_spin.setValue(float(getattr(c, 'tick_minor_width', 0.5)))
            idx_dir = self.tick_direction_combo.findData(
                str(getattr(c, 'tick_direction', 'in') or 'in'))
            if idx_dir >= 0:
                self.tick_direction_combo.setCurrentIndex(idx_dir)
            self.legend_check.setChecked(bool(c.show_legend))
            self.legend_loc_combo.setCurrentText(c.legend_loc)
            self.swatch_check.setChecked(bool(c.show_combo_swatch))
            self.swatch_loc_combo.setCurrentText(c.combo_swatch_loc)
            self.swatch_labels_check.setChecked(bool(c.combo_swatch_labels))
            self.swatch_offset_spin.setValue(float(getattr(c, 'combo_swatch_label_offset', 0.62)))
            self.swatch_inset_scale_spin.setValue(float(getattr(c, 'combo_swatch_inset_scale', 0.24)))
            self.swatch_size_spin.setValue(int(getattr(c, 'combo_swatch_size', 320) or 320))
            self.band_labels_check.setChecked(bool(c.show_band_labels))
            self.band_labels_loc_combo.setCurrentText(c.band_labels_loc)
            self.compass_check.setChecked(bool(c.show_compass))
            self.compass_loc_combo.setCurrentText(c.compass_loc)
            self.beam_check.setChecked(bool(c.show_beam))
            self.beam_loc_combo.setCurrentText(c.beam_loc)
            self.beam_style_combo.setCurrentText(c.beam_style)
            self.scale_bar_check.setChecked(bool(c.show_scale_bar))
            self.scale_bar_asec_spin.setValue(float(c.scale_bar_asec or 0.0))
            idx = self.scale_bar_loc_combo.findData(int(c.scale_bar_loc))
            if idx >= 0:
                self.scale_bar_loc_combo.setCurrentIndex(idx)
            self.scale_bar_color_edit.setText(str(getattr(c, 'scale_bar_color', 'white')))
            self.scale_bar_stroke_edit.setText(str(getattr(c, 'scale_bar_stroke_color', 'black')))
            self.scale_bar_stroke_lw_spin.setValue(float(getattr(c, 'scale_bar_stroke_lw', 1.75)))
        finally:
            for wdg in widgets:
                wdg.blockSignals(False)
        # Dependent UI state that the (now-suppressed) slots would normally set.
        self.combine_blend_combo.setEnabled(c.combine_mode in ('lab', 'hsv', 'hsl'))
        is_custom = self.combine_bg_combo.currentData() == 'custom'
        self.combine_bg_btn.setVisible(is_custom)
        if is_custom:
            self._set_combine_bg_swatch(c.combine_background)
        self._set_tickcolor_swatch(c.tickcolor)
        self._set_facecolor_swatch(c.facecolor)
        self._set_scale_bar_color_swatch(getattr(c, 'scale_bar_color', 'white'))
        self._set_scale_bar_stroke_swatch(getattr(c, 'scale_bar_stroke_color', 'black'))

    # ---------- compose control tabs ----------

    def _build_compose_tabs(self):
        """Group the combined-image controls into a compact tabbed panel."""
        tabs = QTabWidget()

        # --- Compositing: gamma, inverse, color mixing, blend ---
        w = QWidget(); lay = QHBoxLayout(w)
        lay.addWidget(QLabel('Gamma'))
        self.gamma_spin = QDoubleSpinBox()
        self.gamma_spin.setRange(0.1, 10.)
        self.gamma_spin.setSingleStep(0.1)
        self.gamma_spin.setValue(self.session.compose.gamma)
        self.gamma_spin.setToolTip(
            'Compose gamma for the combined image. Single-panel previews look '
            'the same at any gamma (by design). Try 1.5 vs 3.0 to see RYB/Lab changes.')
        self.gamma_spin.valueChanged.connect(self.on_gamma)
        lay.addWidget(self.gamma_spin)
        self.inverse_check = QCheckBox('Inverse (white background)')
        self.inverse_check.toggled.connect(self.on_inverse)
        lay.addSpacing(12)
        lay.addWidget(self.inverse_check)
        lay.addSpacing(16)
        lay.addWidget(QLabel('Color mixing'))
        self.combine_mode_combo = QComboBox()
        self.combine_mode_combo.setToolTip(
            'How overlapping layer colors combine. RGB is the classic channel sum; '
            'CIE Lab mixes in a perceptual space and avoids washing out to white; '
            'RYB mixes subtractively like paint (yellow+blue=green); CMYK like print '
            'ink. HSV/HSL are experimental.')
        for label, value in [('RGB (classic)', 'rgb'), ('CIE Lab', 'lab'),
                             ('HSV (experimental)', 'hsv'), ('HSL (experimental)', 'hsl'),
                             ('RYB (paint)', 'ryb'), ('CMYK (print)', 'cmyk')]:
            self.combine_mode_combo.addItem(label, value)
        self.combine_mode_combo.currentIndexChanged.connect(self.on_combine_mode)
        lay.addWidget(self.combine_mode_combo)
        lay.addWidget(QLabel('Blend'))
        self.combine_blend_combo = QComboBox()
        self.combine_blend_combo.setToolTip(
            'How layer brightnesses combine (Lab/HSV/HSL only). screen accumulates '
            'gracefully; sum is closest to classic RGB; max shows the brightest '
            'layer; mean averages.')
        self.combine_blend_combo.addItems(['screen', 'sum', 'max', 'mean'])
        self.combine_blend_combo.setEnabled(False)
        self.combine_blend_combo.currentTextChanged.connect(self.on_combine_blend)
        lay.addWidget(self.combine_blend_combo)
        lay.addSpacing(16)
        lay.addWidget(QLabel('Background'))
        self.combine_bg_combo = QComboBox()
        self.combine_bg_combo.setToolTip(
            'What shows at zero signal in the combined image. Black is classic; '
            'White is publication-style; Transparent exports a transparent image '
            '(e.g. for slides); or pick any custom color. A non-black background '
            'supersedes the legacy Inverse checkbox.')
        for label, value in [('Black (classic)', 'black'), ('White', 'white'),
                             ('Transparent', 'transparent'), ('Custom color\u2026', 'custom')]:
            self.combine_bg_combo.addItem(label, value)
        self.combine_bg_combo.currentIndexChanged.connect(self.on_combine_bg)
        lay.addWidget(self.combine_bg_combo)
        self.combine_bg_btn = QPushButton()
        self.combine_bg_btn.setFixedWidth(36)
        self.combine_bg_btn.setToolTip('Pick the custom compositing background color')
        self.combine_bg_btn.clicked.connect(self.pick_combine_bg)
        self.combine_bg_btn.setVisible(False)
        lay.addWidget(self.combine_bg_btn)
        lay.addStretch()
        tabs.addTab(w, 'Compositing')

        # --- Colors: palette picker + colorblind check ---
        w = QWidget(); lay = QHBoxLayout(w)
        lay.addWidget(QLabel('Palette'))
        self.palette_combo = QComboBox()
        for group, items in list_palette_menu():
            for label, key in items:
                self.palette_combo.addItem(label, key)
            self.palette_combo.insertSeparator(self.palette_combo.count())
        # trailing separator from last group
        if self.palette_combo.count() and self.palette_combo.itemData(self.palette_combo.count() - 1) is None:
            self.palette_combo.removeItem(self.palette_combo.count() - 1)
        self.palette_combo.setToolTip(
            'Apply a curated or geometry-based palette to the loaded layers, in '
            'order. "Even (auto-N)" generates evenly-spaced colors for however '
            'many layers are loaded.')
        self.palette_combo.currentIndexChanged.connect(self._update_palette_swatch)
        lay.addWidget(self.palette_combo)
        self.palette_swatch = QLabel()
        self.palette_swatch.setFixedHeight(18)
        lay.addWidget(self.palette_swatch)
        self.apply_palette_btn = QPushButton('Apply to layers')
        self.apply_palette_btn.setToolTip("Set each loaded layer's color from this palette")
        self.apply_palette_btn.clicked.connect(self.apply_palette)
        lay.addWidget(self.apply_palette_btn)
        self.tune_colors_btn = QPushButton('Tune colors…')
        self.tune_colors_btn.setToolTip(
            'Open an interactive hue-wheel tuner: spin layer colors while '
            'preserving their spacing, with a live overlapping-circle swatch preview')
        self.tune_colors_btn.clicked.connect(self.open_hue_wheel)
        lay.addWidget(self.tune_colors_btn)
        self.cvd_label = QLabel()
        self.cvd_label.setToolTip('Colour-vision-deficiency check of the current layer colors')
        lay.addWidget(self.cvd_label)
        lay.addStretch()
        tabs.addTab(w, 'Colors')

        # --- Axes & coordinates ---
        w = QWidget()
        outer_axes = QVBoxLayout(w)
        row1 = QHBoxLayout()
        row1.addWidget(QLabel('Coordinates'))
        self.coord_combo = QComboBox()
        self.coord_combo.addItems(['Sexagesimal', 'Decimal'])
        self.coord_combo.currentTextChanged.connect(self.on_coord_style)
        row1.addWidget(self.coord_combo)
        self.minorticks_check = QCheckBox('Minor ticks')
        self.minorticks_check.setChecked(True)
        self.minorticks_check.toggled.connect(self.on_minorticks)
        row1.addWidget(self.minorticks_check)
        row1.addSpacing(16)
        row1.addWidget(QLabel('Tick color'))
        self.tickcolor_btn = QPushButton()
        self.tickcolor_btn.setFixedWidth(36)
        self.tickcolor_btn.clicked.connect(self.pick_tickcolor)
        row1.addWidget(self.tickcolor_btn)
        self.tickcolor_edit = QLineEdit(self.session.compose.tickcolor)
        self.tickcolor_edit.setFixedWidth(80)
        self.tickcolor_edit.editingFinished.connect(self.on_tickcolor_text)
        row1.addWidget(self.tickcolor_edit)
        row1.addStretch()
        outer_axes.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel('X label'))
        self.xlabel_edit = QLineEdit()
        self.xlabel_edit.setPlaceholderText('(auto from FITS)')
        self.xlabel_edit.setMinimumWidth(140)
        self.xlabel_edit.editingFinished.connect(self.on_xlabel)
        row2.addWidget(self.xlabel_edit)
        row2.addSpacing(12)
        row2.addWidget(QLabel('Y label'))
        self.ylabel_edit = QLineEdit()
        self.ylabel_edit.setPlaceholderText('(auto from FITS)')
        self.ylabel_edit.setMinimumWidth(140)
        self.ylabel_edit.editingFinished.connect(self.on_ylabel)
        row2.addWidget(self.ylabel_edit)
        row2.addStretch()
        outer_axes.addLayout(row2)

        row3 = QHBoxLayout()
        for label, attr, default, step in [
            ('Major len', 'tick_major_size_spin', 8.0, 0.5),
            ('Minor len', 'tick_minor_size_spin', 4.0, 0.5),
            ('Major width', 'tick_major_width_spin', 1.0, 0.1),
            ('Minor width', 'tick_minor_width_spin', 0.5, 0.1),
        ]:
            row3.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setRange(0.0, 50.0)
            spin.setSingleStep(step)
            spin.setDecimals(2)
            spin.setValue(default)
            spin.setFixedWidth(64)
            spin.valueChanged.connect(self.on_tick_style)
            setattr(self, attr, spin)
            row3.addWidget(spin)
        row3.addSpacing(12)
        row3.addWidget(QLabel('Tick dir'))
        self.tick_direction_combo = QComboBox()
        for d in ('in', 'out', 'inout'):
            self.tick_direction_combo.addItem(d, d)
        self.tick_direction_combo.setToolTip(
            'Tick direction on the combined WCS plot. WCS axes only support in/out '
            'internally; inout is approximated as inward ticks.')
        self.tick_direction_combo.currentIndexChanged.connect(self.on_tick_direction)
        row3.addWidget(self.tick_direction_combo)
        row3.addStretch()
        outer_axes.addLayout(row3)
        tabs.addTab(w, 'Axes')

        # --- Canvas & annotation: canvas color, title, legend ---
        w = QWidget(); outer = QVBoxLayout(w); lay = QHBoxLayout()
        lay.addWidget(QLabel('Canvas color'))
        self.facecolor_btn = QPushButton()
        self.facecolor_btn.setFixedWidth(36)
        self.facecolor_btn.setToolTip('Background behind the tick labels / axes of the '
                                      'combined plot. Use "none" for a transparent background.')
        self.facecolor_btn.clicked.connect(self.pick_facecolor)
        lay.addWidget(self.facecolor_btn)
        self.facecolor_edit = QLineEdit(self.session.compose.facecolor)
        self.facecolor_edit.setFixedWidth(80)
        self.facecolor_edit.setToolTip("e.g. 'white', 'black', '#204060', or 'none' for transparent")
        self.facecolor_edit.editingFinished.connect(self.on_facecolor_text)
        lay.addWidget(self.facecolor_edit)
        lay.addSpacing(16)
        lay.addWidget(QLabel('Title'))
        self.title_edit = QLineEdit()
        self.title_edit.setMinimumWidth(160)
        self.title_edit.editingFinished.connect(self.on_title)
        lay.addWidget(self.title_edit)
        lay.addSpacing(16)
        self.bare_plot_check = QCheckBox('Bare image (no frame)')
        self.bare_plot_check.setToolTip(
            'Hide ticks, coordinate labels, title, legend, and overlays — image '
            'pixels only (publication / image-only style). Pair with canvas color "none" for '
            'a transparent margin.')
        self.bare_plot_check.toggled.connect(self.on_bare_plot)
        lay.addWidget(self.bare_plot_check)
        lay.addSpacing(16)
        self.legend_check = QCheckBox('Legend')
        self.legend_check.setToolTip('Show a legend mapping each color to its image '
                                     'on the combined figure')
        self.legend_check.toggled.connect(self.on_legend_toggled)
        lay.addWidget(self.legend_check)
        self.legend_loc_combo = QComboBox()
        self.legend_loc_combo.addItems(['upper right', 'upper left',
                                        'lower right', 'lower left'])
        self.legend_loc_combo.setToolTip('Where to place the channel legend')
        self.legend_loc_combo.currentTextChanged.connect(self.on_legend_loc)
        lay.addWidget(self.legend_loc_combo)
        lay.addSpacing(16)
        self.swatch_check = QCheckBox('Color-combo swatch')
        self.swatch_check.setToolTip('Overlapping colored circles showing how the '
                                     'channel colors combine (composite-mode-aware: '
                                     'overlaps use the same RGB/RYB/Lab mixing as the image)')
        self.swatch_check.toggled.connect(self.on_swatch_toggled)
        lay.addWidget(self.swatch_check)
        self.swatch_loc_combo = QComboBox()
        self.swatch_loc_combo.addItems(['lower right', 'lower left',
                                        'upper right', 'upper left'])
        self.swatch_loc_combo.setToolTip('Where to place the color-combination swatch')
        self.swatch_loc_combo.currentTextChanged.connect(self.on_swatch_loc)
        lay.addWidget(self.swatch_loc_combo)
        self.swatch_labels_check = QCheckBox('swatch labels')
        self.swatch_labels_check.setToolTip('Label each circle in the swatch with its channel name')
        self.swatch_labels_check.toggled.connect(self.on_swatch_labels)
        lay.addWidget(self.swatch_labels_check)
        lay.addSpacing(16)
        self.band_labels_check = QCheckBox('Band labels')
        self.band_labels_check.setToolTip("Print each channel's label in its color "
                                          'in a corner of the combined figure')
        self.band_labels_check.toggled.connect(self.on_band_labels_toggled)
        lay.addWidget(self.band_labels_check)
        self.band_labels_loc_combo = QComboBox()
        self.band_labels_loc_combo.addItems(['upper left', 'upper right',
                                             'lower left', 'lower right'])
        self.band_labels_loc_combo.setToolTip('Where to place the colored band labels')
        self.band_labels_loc_combo.currentTextChanged.connect(self.on_band_labels_loc)
        lay.addWidget(self.band_labels_loc_combo)
        lay.addStretch()
        outer.addLayout(lay)

        lay2 = QHBoxLayout()
        lay2.addWidget(QLabel('Swatch label offset'))
        self.swatch_offset_spin = QDoubleSpinBox()
        self.swatch_offset_spin.setRange(0.0, 1.0)
        self.swatch_offset_spin.setSingleStep(0.05)
        self.swatch_offset_spin.setDecimals(2)
        self.swatch_offset_spin.setValue(float(getattr(self.session.compose,
                                                       'combo_swatch_label_offset', 0.62)))
        self.swatch_offset_spin.setToolTip('Push swatch labels outward into each circle outer lobe')
        self.swatch_offset_spin.valueChanged.connect(self.on_swatch_label_offset)
        lay2.addWidget(self.swatch_offset_spin)
        lay2.addSpacing(16)
        lay2.addWidget(QLabel('Swatch scale'))
        self.swatch_inset_scale_spin = QDoubleSpinBox()
        self.swatch_inset_scale_spin.setRange(0.08, 0.45)
        self.swatch_inset_scale_spin.setSingleStep(0.02)
        self.swatch_inset_scale_spin.setDecimals(2)
        self.swatch_inset_scale_spin.setValue(float(getattr(self.session.compose,
                                                            'combo_swatch_inset_scale', 0.24)))
        self.swatch_inset_scale_spin.setToolTip('Corner inset size as a fraction of the plot')
        self.swatch_inset_scale_spin.valueChanged.connect(self.on_swatch_inset_scale)
        lay2.addWidget(self.swatch_inset_scale_spin)
        lay2.addSpacing(8)
        lay2.addWidget(QLabel('Swatch px'))
        self.swatch_size_spin = QSpinBox()
        self.swatch_size_spin.setRange(128, 640)
        self.swatch_size_spin.setSingleStep(32)
        self.swatch_size_spin.setValue(int(getattr(self.session.compose,
                                                   'combo_swatch_size', 320) or 320))
        self.swatch_size_spin.setToolTip('Pixel resolution of the rendered swatch image')
        self.swatch_size_spin.valueChanged.connect(self.on_swatch_size)
        lay2.addWidget(self.swatch_size_spin)
        lay2.addSpacing(16)
        self.compass_check = QCheckBox('Compass')
        self.compass_check.setToolTip('North/east compass rose (requires multicolorfits[overlays] / skyplothelper)')
        self.compass_check.toggled.connect(self.on_compass_toggled)
        lay2.addWidget(self.compass_check)
        self.compass_loc_combo = QComboBox()
        self.compass_loc_combo.addItems(['lower left', 'lower right', 'upper left', 'upper right'])
        self.compass_loc_combo.currentTextChanged.connect(self.on_compass_loc)
        lay2.addWidget(self.compass_loc_combo)
        lay2.addSpacing(8)
        self.beam_check = QCheckBox('Beam')
        self.beam_check.setToolTip('FITS beam ellipse from BMAJ/BMIN (requires [overlays]; header must have beam cards)')
        self.beam_check.toggled.connect(self.on_beam_toggled)
        lay2.addWidget(self.beam_check)
        self.beam_loc_combo = QComboBox()
        self.beam_loc_combo.addItems(['lower left', 'lower right', 'upper left', 'upper right'])
        self.beam_loc_combo.currentTextChanged.connect(self.on_beam_loc)
        lay2.addWidget(self.beam_loc_combo)
        self.beam_style_combo = QComboBox()
        self.beam_style_combo.addItems(['crosshair', 'ellipse', 'crosshairgrid'])
        self.beam_style_combo.currentTextChanged.connect(self.on_beam_style)
        lay2.addWidget(self.beam_style_combo)
        lay2.addSpacing(8)
        self.scale_bar_check = QCheckBox('Scale bar')
        self.scale_bar_check.setToolTip('Arcsecond scale bar (requires multicolorfits[overlays] / skyplothelper)')
        self.scale_bar_check.toggled.connect(self.on_scale_bar_toggled)
        lay2.addWidget(self.scale_bar_check)
        lay2.addWidget(QLabel('asec'))
        self.scale_bar_asec_spin = QDoubleSpinBox()
        self.scale_bar_asec_spin.setRange(0.0, 99999.0)
        self.scale_bar_asec_spin.setDecimals(1)
        self.scale_bar_asec_spin.setSpecialValueText('auto')
        self.scale_bar_asec_spin.setValue(0.0)
        self.scale_bar_asec_spin.setToolTip('Scale-bar length in arcsec (0 = auto)')
        self.scale_bar_asec_spin.valueChanged.connect(self.on_scale_bar_asec)
        lay2.addWidget(self.scale_bar_asec_spin)
        self.scale_bar_loc_combo = QComboBox()
        for label, code in [('lower right', 4), ('lower left', 3), ('upper left', 2), ('upper right', 1)]:
            self.scale_bar_loc_combo.addItem(label, code)
        self.scale_bar_loc_combo.currentIndexChanged.connect(self.on_scale_bar_loc)
        lay2.addWidget(self.scale_bar_loc_combo)
        lay2.addStretch()
        outer.addLayout(lay2)

        lay3 = QHBoxLayout()
        lay3.addWidget(QLabel('Bar color'))
        self.scale_bar_color_btn = QPushButton()
        self.scale_bar_color_btn.setFixedWidth(36)
        self.scale_bar_color_btn.clicked.connect(self.pick_scale_bar_color)
        lay3.addWidget(self.scale_bar_color_btn)
        self.scale_bar_color_edit = QLineEdit(self.session.compose.scale_bar_color)
        self.scale_bar_color_edit.setFixedWidth(72)
        self.scale_bar_color_edit.editingFinished.connect(self.on_scale_bar_color_text)
        lay3.addWidget(self.scale_bar_color_edit)
        lay3.addSpacing(12)
        lay3.addWidget(QLabel('Stroke'))
        self.scale_bar_stroke_btn = QPushButton()
        self.scale_bar_stroke_btn.setFixedWidth(36)
        self.scale_bar_stroke_btn.clicked.connect(self.pick_scale_bar_stroke)
        lay3.addWidget(self.scale_bar_stroke_btn)
        self.scale_bar_stroke_edit = QLineEdit(self.session.compose.scale_bar_stroke_color)
        self.scale_bar_stroke_edit.setFixedWidth(72)
        self.scale_bar_stroke_edit.editingFinished.connect(self.on_scale_bar_stroke_text)
        lay3.addWidget(self.scale_bar_stroke_edit)
        lay3.addWidget(QLabel('width'))
        self.scale_bar_stroke_lw_spin = QDoubleSpinBox()
        self.scale_bar_stroke_lw_spin.setRange(0.0, 10.0)
        self.scale_bar_stroke_lw_spin.setSingleStep(0.25)
        self.scale_bar_stroke_lw_spin.setDecimals(2)
        self.scale_bar_stroke_lw_spin.setValue(float(self.session.compose.scale_bar_stroke_lw))
        self.scale_bar_stroke_lw_spin.valueChanged.connect(self.on_scale_bar_stroke_lw)
        lay3.addWidget(self.scale_bar_stroke_lw_spin)
        lay3.addStretch()
        outer.addLayout(lay3)
        tabs.addTab(w, 'Canvas')

        return tabs

    # ---------- grid alignment ----------

    def refresh_grid_status(self):
        """Show/hide the grid-mismatch warning bar based on loaded panels."""
        report = self.session.grid_report()
        if report['aligned'] or not report['mismatched']:
            self.grid_bar.setVisible(False)
            return
        which = ', '.join(str(i + 1) for i in report['mismatched'])
        plural = len(report['mismatched']) > 1
        ref = report['reference']
        ref_shape = report['shapes'].get(ref)
        shape_txt = (', %d\u00d7%d' % (ref_shape[1], ref_shape[0])) if ref_shape else ''
        msg = ('Panel%s %s %s share the reference grid (panel %d%s). '
               'Use Align layers… to reproject onto a common grid before combining.'
               % ('s' if plural else '', which, "don't" if plural else "doesn't",
                  ref + 1, shape_txt))
        available = reproject_available()
        if not available:
            msg += '  Install the "reproject" package to enable alignment.'
        self.grid_label.setText(msg)
        self.align_btn.setEnabled(available)
        self.grid_bar.setVisible(True)

    def open_align_dialog(self):
        report = self.session.grid_report()
        active = report['active']
        if len(active) < 2:
            return
        dlg = QDialog(self)
        dlg.setWindowTitle('Align layers (reproject)')
        lay = QVBoxLayout(dlg)

        lay.addWidget(QLabel('Reproject all loaded layers onto:'))
        target_combo = QComboBox()
        labels = [('reference', 'Reference panel grid'), ('icrs', 'ICRS (equatorial)'),
                  ('galactic', 'Galactic'), ('fk5', 'FK5 (J2000)'),
                  ('fk4', 'FK4 (B1950)'), ('ecliptic', 'Ecliptic')]
        labels = [(k, v) for k, v in labels if k == 'reference' or k in ALIGN_FRAMES]
        for key, text in labels:
            target_combo.addItem(text, key)
        lay.addWidget(target_combo)

        lay.addWidget(QLabel('Reference panel (defines center / pixel scale):'))
        ref_combo = QComboBox()
        for i in active:
            ref_combo.addItem('Panel %d' % (i + 1), i)
        ref_combo.setCurrentIndex(active.index(report['reference']))
        lay.addWidget(ref_combo)

        north_chk = QCheckBox('North-up in the target frame')
        lay.addWidget(north_chk)
        crop_chk = QCheckBox('Crop to overlap after reproject')
        lay.addWidget(crop_chk)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText('Align')
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        lay.addWidget(buttons)

        if dlg.exec() == QDialog.Accepted:
            self.do_align(target_combo.currentData(), ref_combo.currentData(),
                          north_up=north_chk.isChecked(),
                          crop='overlap' if crop_chk.isChecked() else 'none')

    def do_align(self, target, reference, north_up=False, crop='none'):
        self.show_status('Reprojecting layers\u2026')
        QApplication.processEvents()
        try:
            result = self.session.align_panels(
                target=target, reference=reference, north_up=north_up, crop=crop)
        except Exception as exc:
            QMessageBox.warning(self, 'Alignment failed', str(exc))
            self.show_status('Alignment failed')
            return
        for pw in self.panel_widgets:
            pw.sync_from_state()
            pw.refresh_preview()
        self.refresh_grid_status()
        if self.figure.axes and self.figure.axes[0].images:
            self.plot_combined()
        self.show_status(result.get('message', 'Layers aligned'))

    # ---------- colors & palettes ----------

    def _make_swatch_pixmap(self, colors, box=15, gap=2):
        if not colors:
            return QPixmap()
        pix = QPixmap(len(colors) * (box + gap), box)
        pix.fill(Qt.transparent)
        p = QPainter(pix)
        for k, c in enumerate(colors):
            try:
                p.fillRect(k * (box + gap), 0, box, box, QColor(c))
            except Exception:
                pass
        p.end()
        return pix

    def _update_palette_swatch(self, *_):
        name = self.palette_combo.currentData()
        if not name or name in ('perceptual', 'suggest', 'auto'):
            self.palette_swatch.setPixmap(QPixmap())
            self.palette_swatch.setText(' auto ')
            return
        n = max(2, len(self.session.active_indices()) or 2)
        if name in PALETTES:
            colors = list(PALETTES[name])[:n]
        else:
            try:
                colors = colors_for_hue_pattern(name, n)
            except KeyError:
                colors = []
        self.palette_swatch.setText('')
        self.palette_swatch.setPixmap(self._make_swatch_pixmap(colors))

    def open_hue_wheel(self):
        if len(self.session.active_indices()) < 1:
            self.show_status('Load at least one image first')
            return
        try:
            dlg = HueWheelDialog(self)
        except ValueError as exc:
            self.show_status(str(exc))
            return
        if dlg.exec() == QDialog.Accepted:
            for pw in self.panel_widgets:
                if pw.panel.in_use:
                    pw.sync_from_state()
                    pw.refresh_preview()
            self.refresh_cvd()
            self._replot_if_shown()
            self.show_status('Layer colors updated')

    def apply_palette(self):
        name = self.palette_combo.currentData()
        res = self.session.apply_palette(name)
        for pw in self.panel_widgets:
            pw.sync_from_state()
            pw.refresh_preview()
        self.refresh_cvd()
        if self.figure.axes and self.figure.axes[0].images:
            self._replot_if_shown()
        n = len(res['colors'])
        self.show_status('Applied %s palette to %d layer(s)' % (name, n) if n
                         else 'No layers loaded')

    def refresh_cvd(self):
        report = self.session.colorblind_report()
        if len(report['colors']) < 2:
            self.cvd_label.setText('')
            return
        if report['ok']:
            self.cvd_label.setText('\u2713 colorblind-safe')
            self.cvd_label.setStyleSheet('color: #3fae5a;')
            return
        bits = []
        for kind, res in report['kinds'].items():
            if res['failures']:
                pairs = ', '.join('%d&%d' % (a, b) for a, b, _ in res['failures'])
                bits.append('%s: %s' % (kind[:6], pairs))
        self.cvd_label.setText('\u26a0 confusable (%s)' % '; '.join(bits))
        self.cvd_label.setStyleSheet('color: #e8a04c;')

    # ---------- helpers ----------

    def show_status(self, msg):
        self.statusBar().showMessage(msg, 8000)

    def _draw_placeholder(self):
        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_facecolor('black')
        ax.text(0.5, 0.5,
                'Load FITS images in the panels on the left,\nthen click "Plot Full Resolution".\n\n'
                'Layers must share a common pixel grid.\n'
                'If they don\'t match, use Align layers… or mcf.reproject_image() / mcf.align_stack().',
                color='w', ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        self.canvas.draw_idle()

    def _set_tickcolor_swatch(self, color):
        try:
            from matplotlib.colors import to_hex as mpl_to_hex
            self.tickcolor_btn.setStyleSheet('background-color: %s;' % mpl_to_hex(color))
        except Exception:
            pass

    def _replot_if_shown(self):
        if not self.session.active_panels():
            return
        if self.live_preview_check.isChecked():
            self.schedule_live_preview()
        elif self.figure.axes and self.figure.axes[0].images:
            self.plot_combined()

    def _replot_axes_if_shown(self):
        """Re-render the full WCS plot (ticks, labels, legend) — not live preview."""
        if not self.session.active_panels():
            return
        if self.figure.axes and self.figure.axes[0].images:
            self.plot_combined()

    def schedule_live_preview(self):
        """Debounced fast preview; a no-op unless Fast preview is enabled."""
        if not self.live_preview_check.isChecked():
            return
        if not self.session.active_panels():
            return
        self._live_timer.start(200)

    def _live_preview(self):
        """Fast, downsampled float32 combined preview (plain axes, no WCS)."""
        if not self.session.active_panels():
            return
        try:
            combined = self.session.render_combined(preview=True, max_size=self.LIVE_PREVIEW_MAX)
        except Exception as exc:
            self.show_status('Fast preview failed: %s' % exc)
            return
        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.imshow(np.clip(np.nan_to_num(combined), 0, 1), origin='lower', interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        if self.session.compose.title:
            ax.set_title(self.session.compose.title)
        self.canvas.draw_idle()
        self._wire_cursor_readout()
        self.show_status('Fast preview (fast, downsampled float32 \u2014 click Plot Full Resolution for full WCS)')

    def on_live_preview_toggled(self, checked):
        if not self.session.active_panels():
            return
        if checked:
            if self._combined_plot_ready:
                self.show_status(
                    'Fast preview is pixel-only (no WCS ticks). '
                    'Turn off or click Plot Full Resolution for the full frame.')
            self._live_preview()
        elif self._combined_plot_ready:
            self.plot_combined()
        else:
            self._draw_placeholder()

    def _refresh_wcs_ticklabels(self):
        """Re-assert WCS coordinate tick labels (needed on some Qt backends)."""
        if not self.figure.axes:
            return
        ax = self.figure.axes[0]
        hdr_ast = self.session.common_header.astropy if self.session.common_header else None
        xl, yl = axis_labels_for_compose(hdr_ast, self.session.compose)
        refresh_wcs_ticklabels(ax, self.session.compose, xlabel=xl, ylabel=yl)

    # ---------- slots ----------

    def on_gamma(self, value):
        self.session.compose.gamma = value
        if self.autorefresh_check.isChecked():
            for pw in self.panel_widgets:
                if pw.panel.in_use:
                    pw.refresh_preview()
        self.schedule_live_preview()
        self.show_status('Gamma changed to %.2f' % value)

    def on_panel_preview_hd(self, checked):
        self.session.compose.panel_preview_hd = checked
        if self.autorefresh_check.isChecked():
            for pw in self.panel_widgets:
                if pw.panel.in_use:
                    pw.refresh_preview()

    def on_inverse(self, checked):
        self.session.compose.inverse = checked
        self._replot_if_shown()

    def on_combine_mode(self, *_):
        mode = self.combine_mode_combo.currentData()
        self.session.compose.combine_mode = mode
        self.combine_blend_combo.setEnabled(mode in ('lab', 'hsv', 'hsl'))
        if mode in ('ryb', 'cmyk') and self.combine_bg_combo.currentData() == 'black':
            self.combine_bg_combo.setCurrentIndex(
                self.combine_bg_combo.findData('white'))  # triggers on_combine_bg
        self._replot_if_shown()
        if mode == 'rgb':
            self.show_status('Color mixing: RGB (classic)')
        elif mode == 'ryb':
            self.show_status('Color mixing: RYB (paint / subtractive)')
        elif mode == 'cmyk':
            self.show_status('Color mixing: CMYK (print / subtractive)')
        else:
            self.show_status('Color mixing: %s / %s' % (mode, self.session.compose.combine_blend))

    def on_combine_blend(self, text):
        self.session.compose.combine_blend = text
        self._replot_if_shown()

    def on_combine_bg(self, *_):
        choice = self.combine_bg_combo.currentData()
        self.combine_bg_btn.setVisible(choice == 'custom')
        if choice == 'custom':
            cur = self.session.compose.combine_background
            if str(cur).strip().lower() in ('black', 'white', 'transparent', 'none', ''):
                cur = '#12243a'
            self.session.compose.combine_background = cur
            self._set_combine_bg_swatch(cur)
        else:
            self.session.compose.combine_background = choice
        self._replot_if_shown()

    def pick_combine_bg(self):
        cur = self.session.compose.combine_background
        start = QColor(cur) if QColor.isValidColor(str(cur)) else QColor('#12243a')
        color = QColorDialog.getColor(start, self, 'Compositing background color')
        if color.isValid():
            self.session.compose.combine_background = color.name()
            self._set_combine_bg_swatch(color.name())
            self._replot_if_shown()

    def _set_combine_bg_swatch(self, color):
        try:
            from matplotlib.colors import to_hex as mpl_to_hex
            self.combine_bg_btn.setStyleSheet('background-color: %s;' % mpl_to_hex(color))
        except Exception:
            self.combine_bg_btn.setStyleSheet('')

    def pick_tickcolor(self):
        color = QColorDialog.getColor(QColor('#e5e5e5'), self, 'Tick color')
        if color.isValid():
            self.session.compose.tickcolor = color.name()
            self.tickcolor_edit.setText(color.name())
            self._set_tickcolor_swatch(color.name())
            self._replot_axes_if_shown()

    def on_tickcolor_text(self):
        self.session.compose.tickcolor = self.tickcolor_edit.text().strip()
        self._set_tickcolor_swatch(self.session.compose.tickcolor)
        self._replot_axes_if_shown()

    def pick_facecolor(self):
        cur = self.session.compose.facecolor
        start = QColor(cur) if QColor.isValidColor(str(cur)) else QColor('#ffffff')
        color = QColorDialog.getColor(start, self, 'Canvas (background) color')
        if color.isValid():
            self.session.compose.facecolor = color.name()
            self.facecolor_edit.setText(color.name())
            self._set_facecolor_swatch(color.name())
            self._replot_axes_if_shown()

    def on_facecolor_text(self):
        self.session.compose.facecolor = self.facecolor_edit.text().strip() or 'none'
        self._set_facecolor_swatch(self.session.compose.facecolor)
        self._replot_axes_if_shown()

    def _set_facecolor_swatch(self, color):
        transparent = str(color).strip().lower() in ('none', 'transparent', '')
        if transparent:
            # Checkerboard-ish hint for transparency
            self.facecolor_btn.setStyleSheet(
                'background-color: qlineargradient(x1:0,y1:0,x2:1,y2:1, '
                'stop:0 #bbb, stop:0.5 #fff, stop:1 #bbb);')
            self.facecolor_btn.setText('\u2205')
            return
        self.facecolor_btn.setText('')
        try:
            from matplotlib.colors import to_hex as mpl_to_hex
            self.facecolor_btn.setStyleSheet('background-color: %s;' % mpl_to_hex(color))
        except Exception:
            self.facecolor_btn.setStyleSheet('')

    def on_coord_style(self, text):
        c = self.session.compose
        c.coord_style = text.lower()
        if c.coord_style.startswith('sex'):
            c.x_format, c.y_format = 'hh:mm:ss.ss', 'dd:mm:ss.ss'
        else:
            c.x_format, c.y_format = 'd.dddddd', 'd.dddddd'
        self._replot_axes_if_shown()

    def on_minorticks(self, checked):
        self.session.compose.minorticks = checked
        self._replot_axes_if_shown()

    def on_title(self):
        self.session.compose.title = self.title_edit.text()
        self._replot_axes_if_shown()

    def on_bare_plot(self, checked):
        self.session.compose.bare_plot = checked
        self._replot_axes_if_shown()

    def on_xlabel(self):
        self.session.compose.xlabel = self.xlabel_edit.text().strip()
        self._replot_axes_if_shown()

    def on_ylabel(self):
        self.session.compose.ylabel = self.ylabel_edit.text().strip()
        self._replot_axes_if_shown()

    def on_tick_style(self, *_):
        c = self.session.compose
        c.tick_major_size = float(self.tick_major_size_spin.value())
        c.tick_minor_size = float(self.tick_minor_size_spin.value())
        c.tick_major_width = float(self.tick_major_width_spin.value())
        c.tick_minor_width = float(self.tick_minor_width_spin.value())
        self._replot_axes_if_shown()

    def on_tick_direction(self, _index=None):
        d = self.tick_direction_combo.currentData()
        if d:
            self.session.compose.tick_direction = d
            self._replot_axes_if_shown()

    def on_legend_toggled(self, checked):
        self.session.compose.show_legend = checked
        self._replot_axes_if_shown()

    def on_legend_loc(self, text):
        self.session.compose.legend_loc = text
        if self.session.compose.show_legend:
            self._replot_axes_if_shown()

    def on_swatch_toggled(self, checked):
        self.session.compose.show_combo_swatch = checked
        self._replot_axes_if_shown()

    def on_swatch_loc(self, text):
        self.session.compose.combo_swatch_loc = text
        if self.session.compose.show_combo_swatch:
            self._replot_axes_if_shown()

    def on_swatch_labels(self, checked):
        self.session.compose.combo_swatch_labels = checked
        if self.session.compose.show_combo_swatch:
            self._replot_axes_if_shown()

    def on_swatch_label_offset(self, value):
        self.session.compose.combo_swatch_label_offset = float(value)
        if self.session.compose.show_combo_swatch and self.session.compose.combo_swatch_labels:
            self._replot_axes_if_shown()

    def on_swatch_inset_scale(self, value):
        self.session.compose.combo_swatch_inset_scale = float(value)
        if self.session.compose.show_combo_swatch:
            self._replot_axes_if_shown()

    def on_swatch_size(self, value):
        self.session.compose.combo_swatch_size = int(value)
        if self.session.compose.show_combo_swatch:
            self._replot_axes_if_shown()

    def on_band_labels_toggled(self, checked):
        self.session.compose.show_band_labels = checked
        self._replot_axes_if_shown()

    def on_band_labels_loc(self, text):
        self.session.compose.band_labels_loc = text
        if self.session.compose.show_band_labels:
            self._replot_axes_if_shown()

    def on_compass_toggled(self, checked):
        self.session.compose.show_compass = checked
        self._replot_axes_if_shown()

    def on_compass_loc(self, text):
        self.session.compose.compass_loc = text
        if self.session.compose.show_compass:
            self._replot_axes_if_shown()

    def on_beam_toggled(self, checked):
        self.session.compose.show_beam = checked
        self._replot_axes_if_shown()

    def on_beam_loc(self, text):
        self.session.compose.beam_loc = text
        if self.session.compose.show_beam:
            self._replot_axes_if_shown()

    def on_beam_style(self, text):
        self.session.compose.beam_style = text
        if self.session.compose.show_beam:
            self._replot_axes_if_shown()

    def on_scale_bar_toggled(self, checked):
        self.session.compose.show_scale_bar = checked
        self._replot_axes_if_shown()

    def on_scale_bar_asec(self, value):
        self.session.compose.scale_bar_asec = float(value)
        if self.session.compose.show_scale_bar:
            self._replot_axes_if_shown()

    def on_scale_bar_loc(self, _index):
        code = self.scale_bar_loc_combo.currentData()
        if code is not None:
            self.session.compose.scale_bar_loc = int(code)
            if self.session.compose.show_scale_bar:
                self._replot_axes_if_shown()

    def plot_combined(self):
        if not self.session.active_panels():
            self.show_status('No fits file loaded yet!')
            return
        if self.live_preview_check.isChecked():
            self.live_preview_check.blockSignals(True)
            self.live_preview_check.setChecked(False)
            self.live_preview_check.blockSignals(False)
        try:
            combined = self.session.render_combined()
        except Exception as exc:
            QMessageBox.warning(self, 'Render failed', str(exc))
            return
        setup_combined_axes(self.figure, self.session, combined=combined)
        ax = self.figure.axes[0]
        hdr_ast = self.session.common_header.astropy if self.session.common_header else None
        xl, yl = axis_labels_for_compose(hdr_ast, self.session.compose)
        self.canvas.draw()
        refresh_wcs_ticklabels(ax, self.session.compose, xlabel=xl, ylabel=yl)
        self.canvas.draw_idle()
        self._wire_cursor_readout()
        self._combined_plot_ready = True
        self.show_status('Combined plot updated')

    def reset_session(self):
        """Unload all panels after confirmation."""
        reply = QMessageBox.question(
            self, 'Reset session',
            'Unload all FITS images and restore default compose settings?\n\n'
            'This cannot be undone.',
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if reply != QMessageBox.Yes:
            return
        self.session.reset()
        self.live_preview_check.blockSignals(True)
        self.live_preview_check.setChecked(False)
        self.live_preview_check.blockSignals(False)
        self._combined_plot_ready = False
        self.sync_from_session()
        for pw in self.panel_widgets:
            pw.sync_from_state()
            pw.refresh_preview()
        self.refresh_grid_status()
        self.refresh_cvd()
        self._draw_placeholder()
        self.show_status('Session reset — all panels cleared')

    def _set_scale_bar_color_swatch(self, color):
        try:
            from matplotlib.colors import to_hex as mpl_to_hex
            self.scale_bar_color_btn.setStyleSheet('background-color: %s;' % mpl_to_hex(color))
        except Exception:
            self.scale_bar_color_btn.setStyleSheet('')

    def _set_scale_bar_stroke_swatch(self, color):
        try:
            from matplotlib.colors import to_hex as mpl_to_hex
            self.scale_bar_stroke_btn.setStyleSheet('background-color: %s;' % mpl_to_hex(color))
        except Exception:
            self.scale_bar_stroke_btn.setStyleSheet('')

    def pick_scale_bar_color(self):
        cur = self.session.compose.scale_bar_color
        start = QColor(cur) if QColor.isValidColor(str(cur)) else QColor('#ffffff')
        color = QColorDialog.getColor(start, self, 'Scale bar color')
        if color.isValid():
            self.session.compose.scale_bar_color = color.name()
            self.scale_bar_color_edit.setText(color.name())
            self._set_scale_bar_color_swatch(color.name())
            self._replot_axes_if_shown()

    def on_scale_bar_color_text(self):
        self.session.compose.scale_bar_color = self.scale_bar_color_edit.text().strip() or 'white'
        self._set_scale_bar_color_swatch(self.session.compose.scale_bar_color)
        self._replot_axes_if_shown()

    def pick_scale_bar_stroke(self):
        cur = self.session.compose.scale_bar_stroke_color
        start = QColor(cur) if QColor.isValidColor(str(cur)) else QColor('#000000')
        color = QColorDialog.getColor(start, self, 'Scale bar stroke color')
        if color.isValid():
            self.session.compose.scale_bar_stroke_color = color.name()
            self.scale_bar_stroke_edit.setText(color.name())
            self._set_scale_bar_stroke_swatch(color.name())
            self._replot_axes_if_shown()

    def on_scale_bar_stroke_text(self):
        self.session.compose.scale_bar_stroke_color = self.scale_bar_stroke_edit.text().strip() or 'black'
        self._set_scale_bar_stroke_swatch(self.session.compose.scale_bar_stroke_color)
        self._replot_axes_if_shown()

    def on_scale_bar_stroke_lw(self, value):
        self.session.compose.scale_bar_stroke_lw = float(value)
        if self.session.compose.show_scale_bar:
            self._replot_axes_if_shown()

    def _wire_cursor_readout(self):
        if self._motion_cid is not None:
            self.canvas.mpl_disconnect(self._motion_cid)
        self._motion_cid = self.canvas.mpl_connect('motion_notify_event', self._on_canvas_motion)

    def _on_canvas_motion(self, event):
        if event.xdata is None or event.ydata is None:
            self._cursor_readout.setText('')
            return
        info = self.session.sample_at_pixel(event.xdata, event.ydata)
        if not info['layers']:
            self._cursor_readout.setText('')
            return
        vals = '  '.join('%s: %.5g' % (L['label'], L['value']) for L in info['layers'])
        sky = ''
        if info.get('sky'):
            sky = '  RA=%.5f, Dec=%.5f' % (info['sky']['ra'], info['sky']['dec'])
        self._cursor_readout.setText('x,y=%d,%d   %s%s' % (info['x'], info['y'], vals, sky))

    def save_session(self):
        dlg = PathDialog(
            self, 'Save session',
            start_dir=self.last_file_dir(),
            filters='JSON (*.json);;All files (*)',
            save=True,
            default_name='multicolorfits_session.json',
        )
        if dlg.exec() != QDialog.Accepted:
            return
        path = dlg.path()
        if not path:
            return
        if not path.lower().endswith('.json'):
            path = path + '.json'
        try:
            self.session.save_state(path)
            self.remember_file_dir(path)
            self.show_status('Session saved to %s' % path)
        except Exception as exc:
            QMessageBox.warning(self, 'Save failed', str(exc))

    def load_session(self):
        dlg = PathDialog(
            self, 'Load session',
            start_dir=self.last_file_dir(),
            filters='JSON (*.json);;All files (*)',
            save=False,
        )
        if dlg.exec() != QDialog.Accepted:
            return
        path = dlg.path()
        if not path:
            return
        try:
            warnings = self.session.load_state(path)
            self.remember_file_dir(path)
            self.sync_from_session()
            self.refresh_grid_status()
            self.refresh_cvd()
            if self.live_preview_check.isChecked():
                self._live_preview()
            else:
                self.plot_combined()
            msg = 'Session loaded from %s' % path
            if warnings:
                msg += '  Warnings: ' + '; '.join(warnings)
            self.show_status(msg)
        except Exception as exc:
            QMessageBox.warning(self, 'Load failed', str(exc))

    def save_image(self):
        if not self.session.active_panels():
            self.show_status('No fits file loaded yet!')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save image', self.last_file_dir(),
            'Images (*.png *.jpg *.pdf *.eps *.svg)')
        if not path:
            return
        try:
            transparent = str(self.session.compose.facecolor).strip().lower() in ('none', 'transparent', '')
            self.figure.savefig(path, dpi=300, bbox_inches='tight',
                                facecolor=self.figure.get_facecolor(), transparent=transparent)
            self.remember_file_dir(path)
            self.show_status('Saved %s' % path)
        except Exception as exc:
            QMessageBox.warning(self, 'Save failed', str(exc))

    def export_transparent_cutout(self):
        if not self.session.active_panels():
            self.show_status('No fits file loaded yet!')
            return
        dlg = QDialog(self)
        dlg.setWindowTitle('Export transparent cutout')
        form = QFormLayout(dlg)
        size_spin = QSpinBox()
        size_spin.setRange(0, 8192)
        size_spin.setSpecialValueText('native')
        size_spin.setValue(400)
        size_spin.setToolTip('Longest side in pixels (0 = full resolution)')
        form.addRow('Max size (px)', size_spin)
        lo_spin = QDoubleSpinBox()
        lo_spin.setRange(0, 100)
        lo_spin.setDecimals(2)
        lo_spin.setValue(55.0)
        form.addRow('Alpha low %', lo_spin)
        hi_spin = QDoubleSpinBox()
        hi_spin.setRange(0, 100)
        hi_spin.setDecimals(2)
        hi_spin.setValue(99.3)
        form.addRow('Alpha high %', hi_spin)
        gamma_spin = QDoubleSpinBox()
        gamma_spin.setRange(0.05, 4.0)
        gamma_spin.setDecimals(2)
        gamma_spin.setSingleStep(0.05)
        gamma_spin.setValue(0.5)
        form.addRow('Alpha gamma', gamma_spin)
        src_combo = QComboBox()
        src_combo.addItem('Luminance', 'luma')
        src_combo.addItem('Max RGB', 'max')
        src_combo.addItem('Mean RGB', 'mean')
        form.addRow('Alpha from', src_combo)
        crop_combo = QComboBox()
        crop_combo.addItem('Auto (tight + pad)', 'auto')
        crop_combo.addItem('Full frame', 'none')
        form.addRow('Crop', crop_combo)
        pad_spin = QDoubleSpinBox()
        pad_spin.setRange(0, 0.5)
        pad_spin.setDecimals(2)
        pad_spin.setSingleStep(0.01)
        pad_spin.setValue(0.05)
        form.addRow('Crop pad', pad_spin)
        matte_combo = QComboBox()
        matte_combo.addItem('None', 'none')
        matte_combo.addItem('Circle', 'circle')
        matte_combo.addItem('Ellipse', 'ellipse')
        form.addRow('Matte', matte_combo)
        soft_spin = QDoubleSpinBox()
        soft_spin.setRange(0, 50)
        soft_spin.setDecimals(1)
        soft_spin.setValue(0)
        form.addRow('Soft edge σ (px)', soft_spin)
        invert_check = QCheckBox('Invert (dark = signal)')
        form.addRow(invert_check)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText('Choose file\u2026')
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        form.addRow(buttons)
        if dlg.exec() != QDialog.Accepted:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save transparent cutout', self.last_file_dir(),
            'PNG (*.png);;TIFF (*.tif *.tiff)')
        if not path:
            return
        size = size_spin.value() or None
        try:
            self.session.export_transparent_cutout(
                size=size,
                savepath=path,
                alpha_lo=lo_spin.value(),
                alpha_hi=hi_spin.value(),
                alpha_gamma=gamma_spin.value(),
                alpha_source=src_combo.currentData(),
                crop=crop_combo.currentData(),
                pad=pad_spin.value(),
                matte=matte_combo.currentData(),
                soft_edge=soft_spin.value(),
                invert=invert_check.isChecked(),
            )
            self.remember_file_dir(path)
            self.show_status('Saved transparent cutout %s' % path)
        except Exception as exc:
            QMessageBox.warning(self, 'Export failed', str(exc))

    def save_fits(self):
        if not self.session.active_panels():
            self.show_status('No fits file loaded yet!')
            return
        targets = self.session.fits_save_targets()
        dlg = QDialog(self)
        dlg.setWindowTitle('Save FITS')
        lay = QVBoxLayout(dlg)
        lay.addWidget(QLabel('What to save:'))
        combo = QComboBox()
        for t in targets:
            combo.addItem(t['label'], (t['kind'], t['index']))
        lay.addWidget(combo)
        hint = QLabel("Combined = display-ready RGB cube.  '2D data' = the (possibly\n"
                      "reprojected) single frame.  'colorized RGB' = one layer's color cube.")
        hint.setWordWrap(True)
        lay.addWidget(hint)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText('Choose file\u2026')
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        lay.addWidget(buttons)
        if dlg.exec() != QDialog.Accepted:
            return
        kind, index = combo.currentData()
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save FITS', self.last_file_dir(), 'FITS (*.fits)')
        if not path:
            return
        try:
            self.session.save_fits_target(kind, index=index, savepath=path)
            self.remember_file_dir(path)
            self.show_status('Saved %s' % path)
        except Exception as exc:
            QMessageBox.warning(self, 'Save failed', str(exc))

    def show_params(self):
        text = self.session.params_text()
        print('\n\n' + text + '\n')  # keep the classic print-to-terminal behavior
        TextViewDialog(self, 'Current plot parameters', text).exec()

    def export_script(self):
        TextViewDialog(self, 'Standalone script for the current state',
                       self.session.export_script(), save_suffix='.py').exec()


def main(session=None):
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow(session=session)
    win.show()
    rc = app.exec()
    # Avoid SystemExit in IPython/Jupyter when the window is closed with the X button.
    if __name__ == '__main__':
        sys.exit(rc)
    return rc


if __name__ == '__main__':
    main()
