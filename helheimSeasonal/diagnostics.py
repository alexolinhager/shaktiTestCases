# Video, stills and time series of each Helheim scenario (runme.py), written in
# ./diagnostics. See ../shaktiDiagnostics.py
import os
import sys
import numpy as np
from read_netCDF import read_netCDF
sys.path.insert(0, '..')
from shaktiDiagnostics import shaktiDiagnostics

scenarios = ['seasonal', 'seasonal_firn', 'enhanced_melt']

# Mesh resolution of the runs ('low' or 'high'), as in runme.py
if 'resolution' not in globals():
    resolution = 'low'

for scenario in scenarios:
    name = 'Helheim_' + scenario + '_' + resolution
    if not os.path.isfile(name + '.nc'):
        print('   ' + name + '.nc not found, skipped')
        continue
    md = read_netCDF(name + '.nc')

    # Effective pressure is also averaged below 900 m elevation, where the
    # seasonal runoff reaches the bed in all scenarios
    region = np.where(md.geometry.surface <= 900)[0]

    shaktiDiagnostics(md, name=name, region=region, regionname='Below 900 m elevation')
