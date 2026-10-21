#!/usr/bin/env python3
"""Check upstream optional-grid semantics with absent unused filenames.

Run with Python 3 + numpy. Fixtures and temporary files stay in this folder.
"""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


screen = load('public_screen', HERE/'screen.py')
regression = load('public_regression', HERE.parent/'regression/run_regression.py')


class OptionalGridTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='screen_test_', dir=HERE)
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)/'input'
        regression.write_inputs(self.work, regression.cases()[0])
        self.input = self.work/'tr_in.txt'
        self.expected = screen.screen(screen.read_input(self.input))

    def omit(self, names):
        for name in names:
            (self.work/(name+'.asc')).unlink()

    def constants(self, values):
        lines = self.input.read_text().splitlines()
        lines[7] = ','.join(map(str, values))+',0,90'
        self.input.write_text('\n'.join(lines)+'\n')

    def check_screen(self):
        self.assertEqual(screen.screen(screen.read_input(self.input)), self.expected)

    def test_single_zone_grid_unused(self):
        self.omit(['zones'])
        self.check_screen()

    def test_constant_soil_depth_grid_unused(self):
        self.constants([2,-1,-1])
        self.omit(['zmax'])
        self.check_screen()

    def test_constant_water_table_grid_unused(self):
        self.constants([-1,2,-1])
        self.omit(['depthwt'])
        self.check_screen()

    def test_zero_background_flux_grid_unused(self):
        self.constants([-1,-1,0])
        self.omit(['rizero'])
        self.check_screen()

    def test_all_optional_grids_unused(self):
        self.constants([2,2,0])
        self.omit(['zones','zmax','depthwt','rizero'])
        self.check_screen()

    def test_constant_rainfall_grids_unused(self):
        lines = self.input.read_text().splitlines()
        i = lines.index('cri')
        lines[i+1] = ','.join(['5e-7']*80)
        self.input.write_text('\n'.join(lines)+'\n')
        for path in self.work.glob('ri[0-9]*.asc'):
            path.unlink()
        self.check_screen()

    def test_required_property_grid_missing_is_error(self):
        self.omit(['zmax'])
        with self.assertRaises(FileNotFoundError):
            screen.read_input(self.input)

    def test_multiple_zones_require_zone_grid(self):
        lines = self.input.read_text().splitlines()
        lines[3] = '1,30,100,2'
        lines[11:11] = ['zone,2', lines[9], lines[10]]
        self.input.write_text('\n'.join(lines)+'\n')
        self.omit(['zones'])
        with self.assertRaises(FileNotFoundError):
            screen.read_input(self.input)


if __name__ == '__main__':
    unittest.main(verbosity=2)
