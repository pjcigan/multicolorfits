"""
Interactive hue-wheel tuner for layer colors (Qt).

Shows the composite-mode color-combination swatch for the loaded layers and
lets the user rotate the palette around the hue wheel while preserving
relative spacing.  Panel colors are applied when the dialog is accepted
(default), so dragging the slider stays responsive.
"""

import math

import numpy as np

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QImage, QPainter, QPen, QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QSlider, QVBoxLayout, QWidget,
)

from ..palettes import rotate_colors_from_base


def _rgba_to_qpixmap(rgba):
    """RGBA float array (origin lower-left) -> QPixmap."""
    rgb = np.clip(rgba[..., :3], 0, 1)[::-1]
    alpha = np.clip(rgba[..., 3], 0, 1)[::-1]
    h, w = rgb.shape[:2]
    arr8 = (rgb * 255).astype(np.uint8)
    a8 = (alpha * 255).astype(np.uint8)
    rgba8 = np.dstack([arr8, a8])
    rgba8 = np.ascontiguousarray(rgba8)
    img = QImage(rgba8.data, w, h, 4 * w, QImage.Format_RGBA8888)
    return QPixmap.fromImage(img.copy())


class HueRingWidget(QWidget):
    """Decorative hue ring with a marker for the current rotation."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedSize(120, 120)
        self._rotation = 0.

    def set_rotation(self, deg):
        self._rotation = float(deg) % 360.
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        rect = self.rect().adjusted(8, 8, -8, -8)
        for k in range(36):
            hue = k * 10.
            col = QColor.fromHsv(int(hue * 255 / 360) % 256, 200, 230)
            p.setPen(QPen(col, 6))
            p.drawArc(rect, int((hue - 90) * 16), int(12 * 16))
        cx, cy = rect.center().x(), rect.center().y()
        rad = rect.width() / 2.0
        ang = math.radians(self._rotation - 90)
        mx = cx + rad * math.cos(ang)
        my = cy + rad * math.sin(ang)
        p.setPen(QPen(QColor('#ffffff'), 3))
        p.drawLine(int(cx), int(cy), int(mx), int(my))
        p.end()


class HueWheelDialog(QDialog):
    """
    Rotate active layer hues interactively; preview uses the color-combination swatch.

    By default only the swatch updates while dragging; layer colors are
    committed when OK is clicked.
    """

    def __init__(self, main):
        super().__init__(main)
        self.main = main
        self.session = main.session
        self.base_colors = list(self.session.panel_colors())
        self._labels = [
            self.session.panels[i].label or ('Image %d' % (i + 1))
            for i in self.session.active_indices()
        ]
        if len(self.base_colors) < 1:
            raise ValueError('No loaded layers')

        self._preview_timer = QTimer(self)
        self._preview_timer.setSingleShot(True)
        self._preview_timer.setInterval(50)
        self._preview_timer.timeout.connect(self._refresh_swatch)

        self.setWindowTitle('Tune layer colors')
        self.resize(420, 500)
        lay = QVBoxLayout(self)

        hint = QLabel(
            'Drag the slider to rotate hues around the wheel while keeping their '
            'relative spacing. The overlapping-circle preview uses the same '
            'composite-mode mixing as your combined image. Layer colors are '
            'applied when you click OK.'
        )
        hint.setWordWrap(True)
        lay.addWidget(hint)

        self.swatch_label = QLabel()
        self.swatch_label.setAlignment(Qt.AlignCenter)
        self.swatch_label.setMinimumSize(280, 280)
        self.swatch_label.setStyleSheet('background: #1a1a1a;')
        lay.addWidget(self.swatch_label, stretch=1)

        row = QHBoxLayout()
        self.ring = HueRingWidget()
        row.addWidget(self.ring)
        slider_col = QVBoxLayout()
        slider_col.addWidget(QLabel('Hue rotation'))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 360)
        self.slider.setValue(0)
        self.slider.setToolTip('Rotate the palette around the hue wheel (degrees)')
        self.slider.valueChanged.connect(self._on_slider)
        slider_col.addWidget(self.slider)
        self.deg_label = QLabel('0°')
        self.deg_label.setAlignment(Qt.AlignCenter)
        slider_col.addWidget(self.deg_label)
        row.addLayout(slider_col, stretch=1)
        lay.addLayout(row)

        self.live_panels_check = QCheckBox('Update panel previews while dragging')
        self.live_panels_check.setToolTip(
            'When checked, each loaded panel thumbnail refreshes as you move the '
            'slider (slower on large images). Default updates only the swatch preview.')
        lay.addWidget(self.live_panels_check)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

        self._refresh_swatch()

    def _preview_colors(self, deg=None):
        if deg is None:
            deg = self.slider.value()
        return rotate_colors_from_base(self.base_colors, deg)

    def _on_slider(self, deg):
        self.deg_label.setText('%d°' % deg)
        self.ring.set_rotation(deg)
        if self.live_panels_check.isChecked():
            self._apply_to_panels(deg, refresh_panels=True)
        self._preview_timer.start()

    def _apply_to_panels(self, deg, refresh_panels=False):
        colors = self._preview_colors(deg)
        idx = self.session.active_indices()
        for gi, c in zip(idx, colors):
            self.session.panels[gi].color = c
        if refresh_panels:
            for pw in self.main.panel_widgets:
                if pw.panel.in_use:
                    pw._updating = True
                    try:
                        pw.color_edit.setText(pw.panel.color)
                        pw._set_color_swatch(pw.panel.color)
                    finally:
                        pw._updating = False
                    if self.main.autorefresh_check.isChecked():
                        pw.refresh_preview()
            self.main.refresh_cvd()
            self.main.schedule_live_preview()

    def _refresh_swatch(self):
        data = self.session.render_combo_swatch_colors(
            self._preview_colors(self.slider.value()), labels=self._labels, size=200)
        if data is None:
            self.swatch_label.setPixmap(QPixmap())
            return
        pix = _rgba_to_qpixmap(data['rgba'])
        self.swatch_label.setPixmap(pix.scaled(
            self.swatch_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def accept(self):
        self._apply_to_panels(self.slider.value(), refresh_panels=False)
        super().accept()

    def reject(self):
        idx = self.session.active_indices()
        for gi, c in zip(idx, self.base_colors):
            self.session.panels[gi].color = c
        for pw in self.main.panel_widgets:
            if pw.panel.in_use:
                pw.sync_from_state()
                if self.main.autorefresh_check.isChecked():
                    pw.refresh_preview()
        self.main.refresh_cvd()
        super().reject()
