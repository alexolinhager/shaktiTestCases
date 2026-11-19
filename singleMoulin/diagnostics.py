# Video, stills and time series of the single moulin run (runme.py), written in
# ./diagnostics. See ../shaktiDiagnostics.py
import sys
import numpy as np
from read_netCDF import read_netCDF
sys.path.insert(0, '..')
from shaktiDiagnostics import shaktiDiagnostics

md = read_netCDF('MoulinTransient.nc')

# The moulin, at x=500, y=500
pos = np.argmin(np.sqrt((md.mesh.x - 500)**2 + (md.mesh.y - 500)**2))

# Effective pressure is also averaged downstream of the moulin, between it and
# the outflow at the left edge of the domain
region = np.where(md.mesh.x <= md.mesh.x[pos])[0]

shaktiDiagnostics(md, name='singleMoulin', sites=[[md.mesh.x[pos], md.mesh.y[pos]]],
                  region=region, regionname='Downstream of the moulin')
