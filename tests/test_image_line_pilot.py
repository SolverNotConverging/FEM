"""Geometry and interpolation safeguards for the image-guide pilot."""
import numpy as np
import pytest

from scripts.image_line_pilot import (
    Cell, Config, candidate_valid, channel_flags, contour_segments, merge_intervals,
    single_crossings, sweep_cells,
)


def test_two_strip_internal_and_periodic_gaps():
    intervals = Cell(2, 1.0, 2.0, 3.0).intervals(10.0)
    assert intervals == [(2.0, 3.0), (6.0, 8.0)]
    assert intervals[1][0] - intervals[0][1] == 3.0
    assert intervals[0][0] + 10.0 - intervals[1][1] == 4.0
    assert Cell(1, 2.0).intervals(10.0) == [(4.0, 6.0)]
    with pytest.raises(ValueError, match="outer gap"):
        Cell(2, 4.0, 4.0, 2.0).intervals(10.0)


def test_widths_are_independent_and_single_strip_has_no_internal_gap():
    cells = sweep_cells([1, 2], [1, 3], [.5, 1])
    assert len(cells) == 10
    assert Cell(2, 1, 2, 3) in cells and Cell(2, 2, 1, 3) in cells
    assert all(c.gap_mm == c.w2_mm == 0 for c in cells if c.strips == 1)


def test_contours_recover_known_beta_and_leakage_without_beta_binning():
    widths = [1., 2., 3.]
    x, y = np.meshgrid(widths, widths)
    beta = .1*x + .2*y
    alpha = .01*x + .04*y
    segments = contour_segments(widths, beta, alpha, .6)
    assert segments
    for segment in segments:
        w1, w2, a = segment.T
        np.testing.assert_allclose(.1*w1 + .2*w2, .6)
        np.testing.assert_allclose(a, .01*w1 + .04*w2)
    values = np.concatenate(segments)[:, 2]
    assert values.min() == pytest.approx(.09)
    assert values.max() == pytest.approx(.11)


def test_contours_do_not_bridge_rejected_samples():
    beta = np.array([[.1, .2], [.3, .4]])
    alpha = np.ones((2, 2))*.01
    beta[0, 0] = np.nan
    assert contour_segments([1., 2.], beta, alpha, .25) == []


def test_single_crossings_do_not_bridge_missing_width():
    rows = [dict(strips=1, w1_mm=w, status="accepted", alpha_over_k0=.01*w,
                 beta_harmonic_over_k0=.1*w) for w in (1., 3.)]
    assert single_crossings(rows, [1., 2., 3.], .2) == []
    np.testing.assert_allclose(single_crossings(rows, [1., 3.], .2), [(2., .02)])


def test_disconnected_leakage_ranges_remain_disconnected():
    assert merge_intervals([(3, 4), (0, 1), (.5, 2)]) == [[0, 2], [3, 4]]


def test_backward_growing_and_pml_modes_are_excluded():
    data = dict(direction="forward", alpha_over_k0=.01, pml_fraction=.1, core_fraction=.8)
    assert candidate_valid(data, Config())
    for change in (dict(direction="backward"), dict(alpha_over_k0=-.01),
                   dict(pml_fraction=.8), dict(core_fraction=.01)):
        assert not candidate_valid(data | change, Config())
    signed_seed = data | dict(alpha_over_k0=-2e-5)
    assert candidate_valid(signed_seed, Config(), bound_reference=True)
    assert not candidate_valid(signed_seed, Config())
    assert not candidate_valid(data, Config(), bound_reference=True)


def test_equal_strips_at_half_period_cannot_radiate_odd_harmonic():
    flags = channel_flags(Config(), Cell(2, 1.5, 1.5, 3.5))
    assert flags["half_period_symmetry"] and flags["selected_harmonic_forbidden"]
    assert not channel_flags(Config(), Cell(2, 1., 2., 3.5))["half_period_symmetry"]
    assert not channel_flags(Config(), Cell(2, 1.5, 1.5, 3.0))["half_period_symmetry"]
    assert not channel_flags(Config(harmonic=-2), Cell(2, 1.5, 1.5, 3.5))["selected_harmonic_forbidden"]
