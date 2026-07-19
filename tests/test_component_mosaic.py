"""Tests for component mosaic layout packing and figure construction."""
import numpy as np
import pytest

import multicolorfits as mcf
from multicolorfits.figures import component_slot_grid, make_component_mosaic
from conftest import make_test_data, make_test_header


@pytest.fixture
def multi_session():
    s = mcf.McfSession()
    for i, col in enumerate(['#FF0000', '#00FF00', '#0000FF', '#FFFF00']):
        s.panels[i].set_data(make_test_data(seed=i), make_test_header())
        s.panels[i].color = col
        s.panels[i].label = 'band%d' % i
    return s


class TestComponentSlotGrid:
    def test_single_line(self):
        assert component_slot_grid(3, max_per_line=3, side='top') == [[0, 1, 2]]

    def test_four_with_max_three_top_outer_first(self):
        # Closest to hero fills first → figure order (top): outer then inner.
        grid = component_slot_grid(4, max_per_line=3, side='top')
        assert grid == [
            [3, None, None],  # outer (furthest from hero)
            [0, 1, 2],        # closest to hero
        ]

    def test_four_with_max_three_bottom_closest_first(self):
        grid = component_slot_grid(4, max_per_line=3, side='bottom')
        assert grid == [
            [0, 1, 2],        # closest (just under hero)
            [3, None, None],  # outer
        ]

    def test_max_two_makes_2x2(self):
        grid = component_slot_grid(4, max_per_line=2, side='top')
        assert grid == [[2, 3], [0, 1]]

    def test_left_mirrors_top_packing(self):
        assert component_slot_grid(4, 3, 'left') == component_slot_grid(4, 3, 'top')

    def test_right_mirrors_bottom_packing(self):
        assert component_slot_grid(4, 3, 'right') == component_slot_grid(4, 3, 'bottom')

    def test_empty(self):
        assert component_slot_grid(0, 3, 'top') == []

    def test_bad_side(self):
        with pytest.raises(ValueError):
            component_slot_grid(2, 3, side='diagonal')


class TestMakeComponentMosaic:
    def test_returns_fig_and_axes(self, multi_session):
        fig, axes = make_component_mosaic(multi_session, max_per_line=3, ticks='plain')
        assert fig is not None
        assert axes['combined'] is not None
        assert len(axes['components']) == 4
        assert all(ax is not None for ax in axes['components'])
        # 4 panels, max 3 → 2 strip lines; one empty slot → 1 invisible axes + 4 comps + hero
        assert len(axes['slots']) == 2

    def test_max_per_line_two(self, multi_session):
        fig, axes = make_component_mosaic(
            multi_session, components='top', max_per_line=2, ticks='plain')
        assert len(axes['slots']) == 2
        assert all(len(row) == 2 for row in axes['slots'])

    def test_sides(self, multi_session):
        for side in ('top', 'bottom', 'left', 'right'):
            fig, axes = make_component_mosaic(
                multi_session, components=side, max_per_line=3, ticks='plain')
            assert axes['combined'] is not None
            assert sum(1 for a in axes['components'] if a is not None) == 4

    def test_show_subset_and_order(self, multi_session):
        fig, axes = make_component_mosaic(
            multi_session, show=[2, 0], max_per_line=2, ticks='plain')
        assert len(axes['components']) == 2

    def test_ticks_modes(self, multi_session):
        for mode in ('plain', 'minimal', 'full'):
            fig, axes = make_component_mosaic(
                multi_session, ticks=mode, max_per_line=3)
            assert axes['combined'] is not None

    def test_shared_wcs_projection(self, multi_session):
        fig, axes = make_component_mosaic(multi_session, ticks='minimal')
        hero = axes['combined']
        if hasattr(hero, 'wcs'):
            for ax in axes['components']:
                if ax is not None and hasattr(ax, 'wcs'):
                    assert ax.wcs is not None

    def test_session_wrapper(self, multi_session):
        fig, axes = multi_session.plot_component_mosaic(max_per_line=3)
        assert fig is not None
        assert callable(mcf.make_component_mosaic)

    def test_interior_label_loc(self, multi_session):
        fig, axes = make_component_mosaic(
            multi_session, label_loc='upper left', max_per_line=3)
        ax0 = axes['components'][0]
        assert any(getattr(t, 'get_text', lambda: '')() for t in ax0.texts)

    def test_imagegrid_locator_positions_distinct(self, multi_session):
        """ImageGrid pads via locators — get_position() alone looks collapsed."""
        fig, axes = make_component_mosaic(multi_session, max_per_line=3)
        fig.canvas.draw()
        bboxes = []
        for ax in axes['components']:
            loc = ax.get_axes_locator()
            assert loc is not None  # strip panels come from ImageGrid
            bb = loc(ax, None)
            bboxes.append((round(bb.x0, 4), round(bb.y0, 4),
                           round(bb.width, 4), round(bb.height, 4)))
        assert len(set(bboxes)) == 4

    def test_label_none(self, multi_session):
        fig, axes = make_component_mosaic(
            multi_session, label_loc='none', max_per_line=3)
        ax0 = axes['components'][0]
        assert all(not (getattr(t, 'get_text', lambda: '')() or '').strip()
                   for t in ax0.texts)

    def test_empty_session_raises(self):
        with pytest.raises(ValueError):
            make_component_mosaic(mcf.McfSession())
