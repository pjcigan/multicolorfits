"""suggest_levels picks a starting stretch from the pixel distribution."""

import numpy as np
import pytest

import multicolorfits as mcf


def test_narrow_positive_is_linear():
    rng = np.random.default_rng(0)
    data = 10.0 + rng.normal(0, 0.4, size=(80, 80))
    rec = mcf.suggest_levels(data)
    assert rec['stretch'] == 'linear'
    assert rec['vmin'] < rec['vmax']
    assert rec['span_decades'] is not None and rec['span_decades'] < 1.0
    assert rec['crosses_zero'] is False
    assert 'linear' in rec['reason']


def test_wide_positive_with_hot_tail_clips_vmax():
    # Log-uniform over ~3 decades, plus a thin hot tail above p99.5.
    base = np.logspace(0, 3, 5000)
    data = np.concatenate([base, np.full(20, 1.0e6)])
    rec = mcf.suggest_levels(data, ignore_zeros=False)
    assert rec['stretch'] == 'asinh'
    assert rec['bright_tail'] is True
    assert rec['vmax'] < data.max()
    assert 'saturat' in rec['reason']
    assert rec['span_decades'] > 2.0


def test_few_percent_negative_uses_asinh():
    rng = np.random.default_rng(1)
    data = rng.normal(5.0, 0.5, size=2000)
    data[:100] = -2.0
    rec = mcf.suggest_levels(data)
    assert rec['stretch'] == 'asinh'
    assert rec['crosses_zero'] is True
    assert rec['vmin'] < 0
    assert 'negative' in rec['reason']
    assert 'symlog' in rec['reason']


def test_zeros_change_the_recommendation():
    sky = np.full(400, 10.0)
    gappy = np.concatenate([np.zeros(2000), sky])
    kept = mcf.suggest_levels(gappy, ignore_zeros=False)
    dropped = mcf.suggest_levels(gappy, ignore_zeros=True)
    assert kept['on_floor'] is True
    assert kept['stretch'] == 'asinh'
    assert dropped['stretch'] == 'linear'
    assert dropped['zeros_excluded'] == 2000
    assert 'zeros excluded' in dropped['reason']


def test_constant_does_not_collapse():
    rec = mcf.suggest_levels(np.full((40, 40), 3.0))
    assert rec['stretch'] == 'linear'
    assert rec['vmax'] > rec['vmin']
    assert 'constant' in rec['reason'].lower()


def test_all_nan_raises():
    with pytest.raises(ValueError, match='no finite'):
        mcf.suggest_levels(np.full((8, 8), np.nan))


def test_all_zeros_with_ignore_raises():
    with pytest.raises(ValueError, match='zeros'):
        mcf.suggest_levels(np.zeros(30), ignore_zeros=True)


def test_describe_image_includes_levels(capsys, test_header):
    data = np.full((32, 32), 4.0) + 0.1
    info = mcf.describe_image(data, test_header, name='layer')
    text = capsys.readouterr().out
    assert info['header']['frame']
    assert info['levels']['stretch']
    assert 'Levels suggestion' in text
    assert 'Header' in text


def test_describe_image_verbose_false_is_quiet(capsys, test_header):
    data = np.full((32, 32), 4.0) + 0.1
    info = mcf.describe_image(data, test_header, name='layer', verbose=False)
    assert info['levels']['stretch']
    assert capsys.readouterr().out == ''


def test_describe_images_from_mapping(test_header):
    a = np.linspace(1.0, 2.0, 32 * 32).reshape(32, 32)
    b = np.linspace(10.0, 1000.0, 32 * 32).reshape(32, 32)
    report = mcf.describe_images(
        {'f560w': (a, test_header), 'f770w': (b, test_header)},
        colors=['#00FF00', '#FF0000'],
        verbose=False)
    assert report['names'] == ['f560w', 'f770w']
    assert report['colors'] == ['#00FF00', '#FF0000']
    assert len(report['stretches']) == 2
    assert report['layers']['f560w']['levels']['vmin'] == pytest.approx(report['vmins'][0])
    assert report['layers']['f770w']['color'] == '#FF0000'


def test_describe_images_suggests_colors(test_header):
    data = np.linspace(1.0, 10.0, 32 * 32).reshape(32, 32)
    report = mcf.describe_images([data, data + 1], names=['A', 'B'])
    assert len(report['colors']) == 2
    assert all(c.startswith('#') for c in report['colors'])
    assert report['names'] == ['A', 'B']


def test_panel_apply_suggested_levels():
    s = mcf.McfSession(n_panels=1)
    data = np.linspace(1.0, 2.0, 64 * 64).reshape(64, 64)
    s.panels[0].set_data(data)
    rec = s.panels[0].apply_suggested_levels()
    assert s.panels[0].stretch == rec['stretch']
    assert s.panels[0].vmin == pytest.approx(rec['vmin'])
    assert s.panels[0].vmax == pytest.approx(rec['vmax'])


def test_session_apply_suggested_levels(capsys):
    s = mcf.McfSession(n_panels=2)
    s.panels[0].set_data(np.linspace(1.0, 2.0, 64 * 64).reshape(64, 64))
    s.panels[0].label = 'A'
    s.panels[1].set_data(np.linspace(10.0, 1000.0, 64 * 64).reshape(64, 64))
    recs = s.apply_suggested_levels(verbose=True)
    assert len(recs) == 2
    assert recs[0]['index'] == 0
    assert recs[0]['label'] == 'A'
    assert s.panels[0].stretch == recs[0]['stretch']
    assert s.panels[1].vmin == pytest.approx(recs[1]['vmin'])
    assert 'A:' in capsys.readouterr().out


def test_session_set_display():
    s = mcf.McfSession(n_panels=2)
    data = np.linspace(0.0, 10.0, 32 * 32).reshape(32, 32)
    s.panels[0].set_data(data)
    s.panels[1].set_data(data + 1)
    s.set_display(stretches='asinh', vmins=[0.5, 1.5], vmaxs=8)
    assert s.panels[0].stretch == 'asinh'
    assert s.panels[1].stretch == 'asinh'
    assert s.panels[0].vmin == pytest.approx(0.5)
    assert s.panels[1].vmin == pytest.approx(1.5)
    assert s.panels[1].vmax == pytest.approx(8)
    with pytest.raises(ValueError):
        s.set_display(stretches='nope')
    assert s.panels[0].stretch == 'asinh'
