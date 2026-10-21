#!/usr/bin/env python3
"""Standalone six-case TRIGRS v2.1.0 regression; Python 3 + numpy only.

Finite layer: solve h_t=D h_zz, h_z(0)=-q/Ks, h_z(L)=0 by
propagating its mean and cosine eigenmodes over each constant-flux interval.
Half space: superpose Iverson step responses, evaluated with math.erfc.
Neither reference calls TRIGRS nor imports the project's harness.
Exit 0 = all cases PASS, 1 = numerical FAIL, 2 = invalid input/run/output.
"""
import argparse
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import numpy as np

TOL = 1e-4
L, SLOPE, KS = 2., 35., 1e-5


def cases():
    return [
        ('uniform_past_A5', [.05]*80, .1, False, list(range(81))),
        ('start_burst_then_dry', [1.]*6+[0.]*74, .1, False, list(range(81))),
        ('long_dry_lead_2h_burst', [0.]*80+[.5]*2+[0.]*3, 1., False, list(range(86))),
        ('two_bursts', [.5]*6+[0.]*64+[.5]*6+[0.]*24, .1, False, list(range(101))),
        ('infinite_depth_start_burst', [1.]*6+[0.]*66, .1, True, list(range(73))),
        ('single_output_savage', [1.]*6+[0.]*74, .1, False, [80]),
    ]


def write_inputs(work, case):
    name, rates, F, infinite, times = case
    work.mkdir(parents=True, exist_ok=False)
    theta = SLOPE*math.pi/180
    d0 = F*L**2/3600*math.cos(theta)**2
    def grid(name, value):
        (work/(name+'.asc')).write_text('ncols 1\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 1\nNODATA_value -9999\n'+f'{value:.8g}\n')
    for name, val in [('dem',100),('slope',SLOPE),('zones',1),('zmax',L),('depthwt',L),('rizero',0)]:
        grid(name,val)
    lines = []
    def pair(label, val):
        lines.extend([label,str(val)])
    n = len(rates)
    pair('Name of project', 'Corrected kernel public regression')
    pair('tx,nmax,mmax,zones', f'1,30,{-100 if infinite else 100},1')
    pair('nzs,zmin,uww,nper,t', f'10,0.001,9810,{n},{n*3600}')
    pair('zmax,depth,rizero,min_slope,max_slope','-1,-1,-1,0,90')
    lines.append('zone,1')
    pair('cohesion,phi,uws,diffus,Ksat,Theta_sat,Theta_res,Alpha', f'4000,32,18000,{d0:.17g},{KS},0.45,0.05,-1')
    pair('cri', ','.join(['-1']*n))
    pair('capt', ','.join(str(i*3600) for i in range(n+1)))
    for label,name in [('slope grid','slope'),('elevation grid','dem'),('zones grid','zones'),('depth grid','zmax'),('water table grid','depthwt'),('background flux grid','rizero')]:
        pair(label,name+'.asc')
    lines.append('Rain intensity grids')
    for i,r in enumerate(rates,1):
        grid(f'ri{i:03d}',r*KS)
        lines.append(f'ri{i:03d}.asc')
    for label in ('downstream grid','ordered cells','receptors','weights'):
        pair(label,'no_input')
    pair('Output folder','./')
    pair('Suffix','regression')
    for label in ('Save runoff','Save FS','Save FS depth','Save pressure'):
        pair(label,'F')
    pair('Save water table','F,depth')
    pair('Save infiltration','F')
    pair('Save basal flux','F')
    pair('List detailed output','-2,1')
    pair('Number output times',len(times))
    pair('Output times',','.join(str(t*3600) for t in times))
    pair('Skip other timesteps','T')
    for label,val in [('Analytic fillable porosity','T'),('Positive pressure in rising WT','T'),('psi0=-1/alpha','F'),('Log mass balance','F'),('Flow direction','gener'),('Add steady background to transient','T'),('ASC extension','T'),('Ignore negative pressure for FS','T'),('Ignore capillary fringe','T'),('Deep pressure option','-50,flow')]:
        pair(label,val)
    (work/'tr_in.txt').write_text('\n'.join(lines)+'\n')
    return d0


def reference(case, d0, modes=1024):
    _, rates, _, infinite, outputs = case
    # Use the precision of serialized inputs. References use double arithmetic.
    theta = float(np.float32(np.float32(SLOPE)*np.float32(math.pi/180)))
    phi = float(np.float32(np.float32(32)*np.float32(math.pi/180)))
    D = float(np.float32(d0))/math.cos(theta)**2
    beta = math.cos(theta)**2
    zmin = float(np.float32(.001))
    z = zmin+(L-zmin)*np.arange(11)/10
    flux = np.array([float(f'{r*KS:.8g}') for r in rates],dtype=np.float32)/np.float32(KS)
    flux = np.minimum(flux,1).astype(float)
    n = np.arange(1,modes+1,dtype=float)
    lam = D*(n*math.pi/L)**2
    weights = 2*L*np.cos(z[:,None]*n*math.pi/L)/(n*math.pi)**2
    polynomial = (L-z)**2/(2*L)-L/6
    state = np.zeros(modes)
    mean = 0.
    heads = []
    def iverson(age):
        if age<=0:
            return np.zeros_like(z)
        root = math.sqrt(D*age)
        x = z/(2*root)
        return 2*root*(np.exp(-x*x)/math.sqrt(math.pi)-x*np.array([math.erfc(v) for v in x]))
    for hour in range(len(rates)+1):
        if infinite:
            h = np.zeros_like(z)
            for k in range(min(hour,len(rates))):
                if flux[k]:
                    h += flux[k]*(iverson((hour-k)*3600)-iverson((hour-k-1)*3600))
        elif hour:
            r = flux[hour-1]
            decay = np.exp(-lam*3600)
            state = decay*state+r*(-np.expm1(-lam*3600))
            mean += D*r*3600/L
            h = mean+r*polynomial+weights@(state-r)
        else:
            h = np.zeros_like(z)
        heads.append(np.minimum(beta*(z-L)+h,beta*z))
    p = np.asarray(heads)[outputs]
    denom = 18000*z*math.sin(theta)*math.cos(theta)
    fs = np.minimum(10,np.maximum(math.tan(phi)/math.tan(theta)-9810*np.maximum(p,0)*math.tan(phi)/denom,0)+4000/denom)
    fs[:,z<=.01] = 10
    return z,p,fs


def parse_listing(work, times):
    files = list(work.glob('TRlist_z_p_fs_*.txt'))
    if len(files)!=1:
        raise ValueError('Expected one detailed listing')
    rows, groups, found_times = [], [], []
    for line in files[0].read_text().splitlines():
        if re.match(r'^\s*cell\s+\d+\s',line,re.I) or line.strip().lower().startswith('time:'):
            if rows:
                groups.append(rows)
                rows = []
            found_times.append(float(line.split()[-1].replace('D','E'))/3600)
        elif re.match(r'^\s*[-+.0-9]',line):
            vals = [float(x.replace('D','E')) for x in line.split()]
            if len(vals)==2:  # single-output cell number and angle
                continue
            if len(vals)!=6:
                raise ValueError('Unexpected listing row')
            rows.append(vals)
    if rows:
        groups.append(rows)
    a = np.asarray(groups,dtype=float)
    if a.shape!=(len(times),11,6) or not np.isfinite(a).all():
        raise ValueError(f'Incomplete/nonfinite listing: {a.shape}')
    if found_times and not np.allclose(found_times,times,rtol=0,atol=1e-6):
        raise ValueError('Unexpected listing times')
    return a


def run_case(executable, work, case):
    d0 = write_inputs(work,case)
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['TMPDIR'] = str(work.resolve())
    run = subprocess.run([str(executable)],cwd=work,input='\n',text=True,capture_output=True,env=env,timeout=120)
    # Keep solver diagnostics local; redact the executable's machine path.
    (work/'stdout.log').write_text((run.stdout+run.stderr).replace(str(executable),'<trg>'))
    if run.returncode or 'finished normally' not in (run.stdout+(work/'TrigrsLog.txt').read_text()).lower():
        raise ValueError('TRIGRS did not finish normally')
    a = parse_listing(work,case[-1])
    z,p,fs = reference(case,d0)
    if not np.allclose(a[:,:,0],z,rtol=0,atol=5.1e-5):
        raise ValueError('Unexpected depth nodes')
    # Listing P and FS have five significant figures: rounding is included.
    head_error = float(np.max(np.abs(a[:,:,1]-p)))
    fs_error = float(np.max(np.abs(a[:,:,5]-fs)))
    # Check eigenmode resolution separately, never against the solver.
    if not case[3]:
        _,p2,f2 = reference(case,d0,modes=2048)
        if max(np.max(abs(p-p2)),np.max(abs(fs-f2)))>1e-8:
            raise ValueError('Reference mode resolution insufficient')
    return head_error,fs_error,a


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('executable',type=Path)
    ap.add_argument('--workdir',type=Path,required=True,help='new directory for generated inputs and outputs')
    args = ap.parse_args()
    executable = args.executable.resolve()
    try:
        if not executable.is_file():
            raise ValueError('Executable not found')
        args.workdir.mkdir(parents=True,exist_ok=False)
        print('Independent finite-layer eigenmode recurrence / Iverson half-space')
        print('Tolerance: capped head < 1e-4 m AND FS < 1e-4 (listing rounding included)')
        all_pass = True
        for case in cases():
            he,fe,_ = run_case(executable,args.workdir/case[0],case)
            passed = he<TOL and fe<TOL
            all_pass &= passed
            print(f'{case[0]}: max capped-head error={he:.9g} m; max FS error={fe:.9g}; {"PASS" if passed else "FAIL"}',flush=True)
        print('OVERALL '+('PASS' if all_pass else 'FAIL'))
        return 0 if all_pass else 1
    except (ValueError,OSError,subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':
    sys.exit(main())
