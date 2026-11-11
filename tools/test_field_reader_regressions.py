"""Saved-field reader contracts; synthetic inputs only, no solver execution."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import plot_fields_b5 as p
from field_io import FieldFormatError


def frame(step, rows, extra=""):
    columns = "id xc yc xlo ylo xhi yhi nrho u v temp" + extra
    return (f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF CELLS\n{len(rows)}\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 .1\n0 .02\n-.5 .5\n"
            f"ITEM: CELLS {columns}\n" +
            "".join(" ".join(map(str, row)) + "\n" for row in rows))


def row(n=1, u=10, t=0, ident=1):
    return [ident, .025, .005, .02, 0, .03, .01, n, u, 0, t]


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "field.grid"

    def read(self, text):
        self.path.write_text(text)
        original = self.path.read_bytes()
        result = p.average_tail(self.path, frac=1)
        self.assertEqual(self.path.read_bytes(), original)
        return result

    def test_unequal_window_weights_and_sentinel(self):
        d = self.read(frame(1, [row()]) + frame(2, [row(3, 20, 10)]))
        np.testing.assert_allclose([d['nrho'][0], d['u'][0], d['t'][0]], [2, 17.5, 10])

    def test_both_measured_temperatures(self):
        d = self.read(frame(1, [row(t=2)]) + frame(2, [row(3, 20, 10)]))
        self.assertEqual(d['t'][0], 8)

    def test_no_measured_temperature_retains_missing_sentinel(self):
        d = self.read(frame(1, [row()]) + frame(2, [row(3, 20, 0)]))
        self.assertEqual(d['t'][0], 0)

    def test_missing_temperature_ignored(self):
        d = self.read(frame(1, [row(t=float('nan'))]) + frame(2, [row(3, 20, 10)]))
        self.assertEqual(d['t'][0], 10)

    def test_unsupported_width_rejected(self):
        with self.assertRaises(FieldFormatError):
            self.read(frame(1, [row() + [1, 2]], ' vol other'))

    def test_unequal_flow_area_region_weights(self):
        d = self.read(frame(1, [row(0, 0, 0) + [1], row(4, 20, 10, 2) + [3]], ' vol'))
        self.assertEqual(p.region_mean(d, .01, .055, .017), (3., 2))

    def test_zero_density_window_included(self):
        # Second cell keeps both windows populated, including the empty first cell.
        d = self.read(frame(1, [row(0, 0, 0), row(1, 1, 4, 2)]) +
                      frame(2, [row(2, 20, 10), row(1, 1, 4, 2)]))
        np.testing.assert_allclose([d['nrho'][0], d['u'][0], d['t'][0]], [1, 20, 10])

    def test_flow_area_and_zero_density_region(self):
        rows = [row(0, 0, 0) + [.0001], row(2, 20, 10, 2) + [.0001]]
        d = self.read(frame(1, rows, ' vol'))
        self.assertEqual(p.region_mean(d, .01, .055, .017), (1., 2))

    def test_wrong_extra_column_rejected(self):
        with self.assertRaises(FieldFormatError):
            self.read(frame(1, [row() + [1]], ' other'))

    def test_mixed_layout_rejected(self):
        with self.assertRaises(FieldFormatError):
            self.read(frame(1, [row()]) + frame(2, [row() + [1]], ' vol'))

    def test_invalid_flow_area_rejected(self):
        for area in [-1, float('nan'), float('inf')]:
            with self.subTest(area=area), self.assertRaises(FieldFormatError):
                self.read(frame(1, [row() + [area]], ' vol'))

    def test_changed_flow_area_rejected(self):
        with self.assertRaises(FieldFormatError):
            self.read(frame(1, [row() + [1]], ' vol') + frame(2, [row() + [2]], ' vol'))

    def test_flow_area_alone_does_not_populate_frame(self):
        with self.assertRaises(FieldFormatError):
            self.read(frame(1, [row(0, 0, 0) + [1]], ' vol'))

    def test_reordered_ids_and_incomplete_tail(self):
        d = self.read(frame(1, [row(), row(1, 1, 4, 2)]) +
                      frame(2, [row(1, 1, 4, 2), row(3, 20, 10)]) +
                      'ITEM: TIMESTEP\n3\n')
        self.assertEqual(d['steps'], [1, 2])
        self.assertEqual(d['u'][0], 17.5)


    def test_default_tail_rejects_changed_earlier_grid(self):
        for mutation in ('area', 'negative', 'nan', 'inf', 'id', 'geometry'):
            with self.subTest(mutation=mutation):
                early = [row(0, 0, 0) + [1], row(4, 20, 10, 2) + [3]]
                late = [r.copy() for r in early]
                if mutation == 'area':
                    early[0][-1], early[1][-1] = 3, 1
                elif mutation in ('negative', 'nan', 'inf'):
                    early[0][-1] = {'negative': -1, 'nan': float('nan'), 'inf': float('inf')}[mutation]
                elif mutation == 'id':
                    early[0][0] = 3
                else:
                    early[0][1] += .001
                self.path.write_text(frame(1, early, ' vol') + frame(2, late, ' vol'))
                with self.assertRaises((FieldFormatError, SystemExit)):
                    p.average_tail(self.path)

    def test_default_tail_stable_grid_trend(self):
        early = [row(0, 0, 0) + [1], row(4, 20, 10, 2) + [3]]
        late = [row(4, 20, 10) + [1], row(0, 0, 0, 2) + [3]]
        self.path.write_text(frame(1, early, ' vol') + frame(2, late, ' vol'))
        d = p.average_tail(self.path)
        self.assertEqual(d['steps'], [2])
        self.assertEqual(d['trend'], [(1, 3.), (2, 1.)])

if __name__ == '__main__':
    unittest.main()
