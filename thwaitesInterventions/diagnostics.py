# Video, stills and time series of each Thwaites run (runme.py), written in
# ./diagnostics. See ../shaktiDiagnostics.py
import os
import pickle
import sys
import numpy as np
from scipy.spatial import cKDTree
from read_netCDF import read_netCDF
sys.path.insert(0, '..')
from shaktiDiagnostics import shaktiDiagnostics

runs = ['ThwaitesSpinup', 'ThwaitesCapacity', 'Thwaites_control', 'Thwaites_pump', 'Thwaites_inject']

# Mesh resolution of the runs ('low' or 'high'), as in runme.py
if 'resolution' not in globals():
    resolution = 'low'

# Intervention sites, as in runme.py
pumpXY = np.array([[-1429020, -441873],
                   [-1438900, -424922],
                   [-1419150, -463010]])

for run in runs:
    name = run + '_' + resolution
    if not os.path.isfile(name + '.nc'):
        print('   ' + name + '.nc not found, skipped')
        continue
    md = read_netCDF(name + '.nc')

    # Effective pressure is also averaged within 200 km of the grounding line:
    # the points where md.mask.ocean_levelset is zero along the mesh edges
    el = md.mesh.elements.astype(int) - 1
    edges = np.unique(np.sort(np.vstack((el[:, [0, 1]], el[:, [1, 2]], el[:, [2, 0]])), axis=1), axis=0)
    a = md.mask.ocean_levelset[edges[:, 0]]
    b = md.mask.ocean_levelset[edges[:, 1]]
    cut = np.where((a > 0) != (b > 0))[0]
    f = (a[cut] / (a[cut] - b[cut]))[:, None]
    xy = np.column_stack((md.mesh.x, md.mesh.y))
    gl = xy[edges[cut, 0], :] * (1 - f) + xy[edges[cut, 1], :] * f
    distance = cKDTree(gl).query(xy)[0]
    region = np.where(distance <= 200e3)[0]
    options = {'name': name, 'sites': pumpXY, 'region': region, 'regionname': 'Within 200 km of the grounding line'}

    # Capacity run: the head is pinned over the sites, so what the pumps take
    # out is the net inflow into the sites computed in step 2 of runme.py
    capfile = 'ThwaitesCapacity_cap_' + resolution + '.pkl'
    if run == 'ThwaitesCapacity' and os.path.isfile(capfile):
        with open(capfile, 'rb') as fid:
            cap = pickle.load(fid)
        options['series'] = [{'name': 'Pump output', 'time': cap['time'], 'value': np.sum(cap['Qin'], axis=0)}]

    shaktiDiagnostics(md, **options)
