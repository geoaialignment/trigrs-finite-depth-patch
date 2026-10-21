#!/usr/bin/env python3
"""Input-only S21 and S21-D1 screens for TRIGRS v2.1.0; numpy only.

Grid paths are relative to tr_in.txt, or --base-dir when specified. Reports
zone summaries over actual cell slopes, or both endpoints of --slopes.
Exposure is a branch diagnostic, not a prediction of error or instability.
"""
import argparse
import math
from pathlib import Path
import sys
import numpy as np


def numbers(text):
    return np.asarray([float(x.replace('D','E').replace('d','e')) for x in text.replace(',',' ').split()],dtype=np.float32)


def read_grid(path):
    lines = path.read_text().splitlines()
    header = dict((s.split()[0].lower(),float(s.split()[1])) for s in lines[:6])
    a = np.asarray([float(x) for s in lines[6:] for x in s.split()],dtype=np.float32)
    if a.size!=int(header['ncols']*header['nrows']):
        raise ValueError('Grid size mismatch: '+path.name)
    return a,header.get('nodata_value',-9999)


def read_input(path, base=None):
    lines = [s.strip() for s in path.read_text().splitlines() if s.strip()]
    # Native v2.1.0 reads these fields positionally; printed labels vary.
    tx,_,mmax,nzones = numbers(lines[3])
    settings = numbers(lines[5])
    nper = int(settings[3])
    constants = numbers(lines[7])
    props = {}
    i = 8
    for _ in range(int(nzones)):
        zone = int(numbers(lines[i].split(',',1)[1])[0])
        props[zone] = numbers(lines[i+2])
        i += 3
    cri = numbers(lines[i+1]); bounds = numbers(lines[i+3]).astype(float)
    i += 4
    if len(cri)!=nper or len(bounds)!=nper+1 or (np.diff(bounds)<=0).any():
        raise ValueError('Invalid forcing periods')
    base = base or path.parent
    filenames = {}
    for name in ('slope','dem','zones','L','wt','q0'):
        filenames[name] = base/lines[i+1].strip("'\"")
        i += 2
    grids = {name: read_grid(filenames[name]) for name in ('slope','dem')}
    count = len(grids['slope'][0])
    # Match trigrs.f90: a single zone and nonnegative property constants
    # replace their grids entirely; unused filenames need not exist.
    grids['zones'] = ((np.ones(count,dtype=np.float32),-9999)
                      if nzones==1 else read_grid(filenames['zones']))
    for name,k in [('L',0),('wt',1),('q0',2)]:
        grids[name] = ((np.full(count,constants[k],dtype=np.float32),-9999)
                       if constants[k]>=0 else read_grid(filenames[name]))
    if any(len(g[0])!=count for g in grids.values()):
        raise ValueError('Input grids have different sizes')
    i += 1  # rainfall heading
    rain = []
    for k in range(nper):
        if cri[k]>=0:
            rain.append(np.full(count,cri[k]))
        else:
            r,_ = read_grid(base/lines[i+k].strip("'\""))
            if len(r)!=count:
                raise ValueError('Rainfall grid size mismatch')
            rain.append(r)
    # A scalar rainfall override means the named rainfall grid is unused.
    i += nper
    def locate(test):
        return next(j for j in range(i,len(lines)) if test(lines[j].lower()))
    nout_i = locate(lambda s:s.startswith('number') and ('output' in s or 'save output' in s))
    nout = int(numbers(lines[nout_i+1])[0])
    times = numbers(lines[nout_i+3]).astype(float)
    if len(times)!=nout or (times<0).any():
        raise ValueError('Invalid requested output times')
    valid = np.ones(count,dtype=bool)
    for name in ('slope','dem','zones'):
        a,nd = grids[name]
        valid &= (a!=nd)&np.isfinite(a)
    for name,k in [('L',0),('wt',1),('q0',2)]:
        a,nd = grids[name]
        valid &= (a!=nd)&np.isfinite(a)
    rain = np.asarray(rain).T
    if not np.isfinite(rain[valid]).all() or (rain[valid]<0).any():
        raise ValueError('Missing/negative rainfall in active cells')
    return dict(props=props,grids=grids,valid=valid,rain=rain,bounds=bounds,times=times,finite=mmax>0)


def screen(data, slope_range=None):
    g = data['grids']; report = []
    for zone,properties in sorted(data['props'].items()):
        mask = data['valid']&(g['zones'][0]==zone)
        indices = np.flatnonzero(mask)
        if not len(indices):
            continue
        maxage, exposure, diagnostic, applicable = 0.,False,False,0
        slopes = []
        for j in indices:
            L = float(g['L'][0][j]); wt = float(g['wt'][0][j])
            if L<=0 or properties[3]<=0 or properties[4]<=0:
                raise ValueError('Nonpositive depth, diffusivity or conductivity')
            alpha = float(properties[7])
            app = bool(data['finite'] and (alpha<=0 or wt<=1/alpha))
            applicable += int(app)
            used = slope_range if slope_range is not None else [float(g['slope'][0][j])]
            for slope in used:
                if not 0<=slope<90:
                    raise ValueError('Slope must be in [0,90) degrees')
                slopes.append(slope)
                # Match native single precision storage of angle conversion.
                theta = float(np.float32(np.float32(slope)*np.float32(math.pi/180)))
                D = float(properties[3])/math.cos(theta)**2
                ages = D*(data['times'][:,None]-data['bounds'][None,:-1])/L**2
                maxage = max(maxage,float(np.maximum(ages,0).max()))
                exposure |= bool(app and (ages>5).any())
                if app:
                    wet = np.minimum(data['rain'][j],properties[4])>0
                    for row in ages:
                        for k in np.flatnonzero(row>=72):
                            diagnostic |= bool(np.any((row[k+1:]>0)&wet[k+1:]))
        report.append(dict(zone=zone,cells=len(indices),applicable_cells=applicable,
                           slope_min=min(slopes),slope_max=max(slopes),max_A=maxage,
                           exposed=exposure,dry_lead=diagnostic))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path,help='tr_in.txt')
    p.add_argument('--base-dir',type=Path,help='base for relative grid filenames')
    p.add_argument('--slopes',nargs=2,type=float,metavar=('MIN','MAX'),help='evaluate a slope range for each cell')
    a = p.parse_args()
    try:
        if a.slopes and a.slopes[0]>a.slopes[1]:
            raise ValueError('Reversed slope range')
        rows = screen(read_input(a.input,a.base_dir),a.slopes)
        if not rows:
            raise ValueError('No active property-zone cells')
        print('Requested output times only; period starts include first and zero-rain periods.')
        print('S21: A > 5. S21-D1: A >= 72 AND a younger wet period has started.')
        print('A = (D0 / cos(slope)^2) * (output_time - period_start) / L^2')
        for r in rows:
            label = 'N/A' if not r['applicable_cells'] else ('EXPOSED' if r['exposed'] else 'UNEXPOSED')
            d1 = 'N/A' if not r['applicable_cells'] else ('EXPOSED-D1' if r['dry_lead'] else 'NO-D1')
            print(f"zone {r['zone']}: cells={r['cells']}; applicable={r['applicable_cells']}; slopes={r['slope_min']:.8g}..{r['slope_max']:.8g} deg; max A={r['max_A']:.9g}; S21={label}; S21-D1={d1}")
        print('Applicability: finite saturated kernel or unsaturated fallback d <= 1/alpha.')
        if a.slopes:
            print('Slope range endpoints give age bounds; flags are conservative over the range.')
        return 0
    except (ValueError,OSError,IndexError,StopIteration) as e:
        print('ERROR: '+str(e),file=sys.stderr)
        return 2

if __name__=='__main__':
    sys.exit(main())
