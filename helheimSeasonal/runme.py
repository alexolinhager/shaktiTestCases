# Standalone SHAKTI simulations of Helheim Glacier forced by seasonal meltwater
# inputs to the bed, after the "Seasonal Hydrology Forcing" experiments of
# Sommers et al. (2024, https://doi.org/10.1029/2024GL109168, section 2.4.1).
#
# Unlike in the paper, the stress balance is not solved here: the ice velocity
# stays at its winter value while the subglacial hydrology evolves.
#
# Python version of runme.m
import numpy as np
from scipy.io import loadmat
from generic import generic
from read_netCDF import read_netCDF
from solve import solve
from write_netCDF import write_netCDF

# Meltwater scenarios to run one after the other, each one starts from the
# winter base state and is saved in Helheim_<scenario>_<resolution>.nc:
#    'seasonal'       seasonal runoff below 900 m
#    'seasonal_firn'  seasonal runoff below 900 m + steady firn aquifer drainage
#    'enhanced_melt'  doubled runoff below 1500 m + doubled firn aquifer drainage
scenarios = ['seasonal', 'seasonal_firn', 'enhanced_melt']

# Number of processors and mesh resolution: change them according to your
# machine, or set nproc and resolution before running this script (as
# run/runTestCases.py does)
#    'low'   300 m to 3 km (7,366 elements)
#    'high'  mesh of the paper (27,913 elements)
if 'nproc' not in globals():
    nproc = 2
if 'resolution' not in globals():
    resolution = 'low'

for scenario in scenarios:
    print('Seasonal meltwater forcing and solve (' + scenario + ')')
    if resolution == 'low':
        md = read_netCDF('HelheimWinterSpinup_300m.nc')
    elif resolution == 'high':
        md = read_netCDF('HelheimWinterSpinup_70m.nc')
    else:
        raise ValueError("resolution '" + resolution + "' not supported")

    # The hydrology does not use the enthalpy, but ISSM's Python interface
    # fails to rebuild it from the temperature when the pressure is not set
    md.initialization.enthalpy = md.materials.heatcapacity * (md.initialization.temperature - md.constants.referencetemperature)

    if scenario == 'seasonal':
        melt_factor, melt_elevation, firn_factor = 1, 900, 0
    elif scenario == 'seasonal_firn':
        melt_factor, melt_elevation, firn_factor = 1, 900, 1
    elif scenario == 'enhanced_melt':
        melt_factor, melt_elevation, firn_factor = 2, 1500, 2
    else:
        raise ValueError("scenario '" + scenario + "' not supported")

    # Define the time stepping scheme: run for 1 year with a time step of 6 hr,
    # and save the model output every 10 days
    md.timestepping.start_time = 0
    md.timestepping.time_step = 6 * 3600 / md.constants.yts  # 6 hr Time step (in years)
    md.timestepping.final_time = 1
    md.settings.output_frequency = 4 * 10  # output every 10 days
    time = md.timestepping.time_step * np.arange(round(md.timestepping.final_time / md.timestepping.time_step) + 1)

    # Meltwater input rate (m/yr) from hourly 2018 MERRA-2 runoff near the
    # terminus, starting on January 1st and going back to zero at the end of the year
    runoff = np.append(loadmat('merra2_runoff.mat')['ro2018'].flatten(), 0)
    runoff = melt_factor * np.interp(time, np.arange(runoff.size) / (runoff.size - 1), runoff)

    # Distributed input to the bed (m/yr) at low elevation, where the ice is
    # heavily crevassed
    pos = np.where(md.geometry.surface <= melt_elevation)[0]
    md.hydrology.englacial_input = np.zeros((md.mesh.numberofvertices + 1, time.size))
    md.hydrology.englacial_input[-1, :] = time
    md.hydrology.englacial_input[pos, :] = runoff

    # Steady drainage of the firn aquifer (50e6 m^3/yr) through crevasses at
    # 1500 m elevation, divided evenly between point inputs (m^3/s)
    pos = np.where((md.geometry.surface >= 1500) & (md.geometry.surface <= 1515))[0]
    md.hydrology.moulin_input = np.zeros(md.mesh.numberofvertices)
    md.hydrology.moulin_input[pos] = firn_factor * 50e6 / md.constants.yts / pos.size

    # Specify that you want to run the model on your current computer
    md.cluster = generic('np', nproc)

    md.miscellaneous.name = 'Helheim_' + scenario + '_' + resolution
    md.verbose.solution = True
    md = solve(md, 'Transient')

    write_netCDF(md, 'Helheim_' + scenario + '_' + resolution + '.nc')
