#!/usr/bin/env python3
"""Build serial or MPI TRIGRS in a separate output directory (standard library)."""
import argparse
import os
from pathlib import Path
import subprocess

SOURCES = {'serial': ['grids.f95', 'input_vars.f95', 'model_vars.f95', 'dsimps.f95', 'input_file_defs.f95', 'iverson.f95', 'pstpi.f95', 'satfin.f95', 'savage.f95', 'steady.f95', 'trini.f95', 'unsinf.f95', 'ivestp.f95', 'pstpf.f95', 'rnoff.f95', 'satinf.f95', 'smallt.f95', 'svgstp.f95', 'unsfin.f95', 'unsth.f95', 'ssizgrd.f95', 'svlist.f95', 'trigrs.f90', 'flux.f90', 'prpijz.f90', 'svxmdv.f90', 'svijz.f90', 'dzero_brac.f90', 'calerf.f', 'dbsct.f', 'derfc.f', 'irdgrd.f', 'irdswm.f', 'isvgrd.f', 'roots.f', 'srdgrd.f', 'srdswm.f', 'ssvgrd.f'], 'mpi': ['modules_p.f90', 'grids.f95', 'input_vars.f95', 'model_vars.f95', 'dsimps.f95', 'input_file_defs.f95', 'iverson.f95', 'pstpi.f95', 'pstpi_p.f95', 'satfin.f95', 'satfin_p.f95', 'savage.f95', 'steady.f95', 'trini.f95', 'trini_p.f90', 'unsinf.f95', 'unsinf_p.f95', 'ivestp.f95', 'ivestp_p.f95', 'pstpf.f95', 'pstpf_p.f95', 'rnoff.f95', 'satinf.f95', 'satinf_p.f95', 'smallt.f95', 'svgstp.f95', 'svgstp_p.f95', 'unsfin.f95', 'unsfin_p.f95', 'unsth.f95', 'unsth_p.f95', 'ssizgrd.f95', 'ssizgrd_p.f95', 'svlist.f95', 'rnoff_p.f95', 'steady_p.f95', 'trigrs_p.f90', 'partial_p.f90', 'flux.f90', 'flux_p.f90', 'prpijz.f90', 'svxmdv.f90', 'svijz.f90', 'dzero_brac.f90', 'srdgrd_p.f90', 'irdgrd_p.f90', 'calerf.f', 'dbsct.f', 'derfc.f', 'irdgrd.f', 'irdswm.f', 'irdswm_p.f90', 'isvgrd.f', 'roots.f', 'srdgrd.f', 'srdswm.f', 'srdswm_p.f90', 'ssvgrd.f']}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path, help='patched upstream src directory')
    p.add_argument('output', type=Path, help='build directory')
    p.add_argument('--mode', choices=SOURCES, default='serial')
    p.add_argument('--compiler', help='defaults to gfortran or mpif90 (MPICH)')
    a = p.parse_args()
    source = a.source.resolve()
    a.output.mkdir(parents=True, exist_ok=True)
    compiler = a.compiler or ('gfortran' if a.mode == 'serial' else 'mpif90')
    command = [compiler, '-O2', '-std=legacy', '-fallow-argument-mismatch',
               '-ffree-line-length-none', '-o', 'trg']
    command += [str(source / name) for name in SOURCES[a.mode]]
    env = os.environ.copy()
    env['TMPDIR'] = str(a.output.resolve())
    subprocess.run(command, cwd=a.output, env=env, check=True)
    print('Built', a.mode, 'trg')

if __name__ == '__main__':
    main()
