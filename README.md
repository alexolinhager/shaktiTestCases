# shaktiTestCases

Standard test cases for SHAKTI, the subglacial hydrology model of [ISSM](https://issmteam.github.io/ISSM-Documentation/). Each case runs SHAKTI on its own (the ice does not move or change) and comes in a MATLAB and a Python version.

## Test cases

| Case | What it runs | Source |
|---|---|---|
| `singleMoulin` | One moulin feeding 4 m³/s to the centre of a 1 km square under a 500 m thick tilted slab, for 30 days. A channel grows from the moulin to the outflow edge. | Sommers et al. ([2018](https://doi.org/10.5194/gmd-11-2955-2018)); it is the `shakti` example of ISSM |
| `helheimSeasonal` | Helheim Glacier, Greenland, for 1 year from a winter state, with three meltwater scenarios: seasonal runoff, seasonal runoff plus firn aquifer drainage, and enhanced melt. | "Seasonal Hydrology Forcing" experiments of Sommers et al. ([2024](https://doi.org/10.1029/2024GL109168), section 2.4.1) |
| `thwaitesInterventions` | Thwaites Glacier, Antarctica: a 5-year spin-up, a 5-year run that measures how much water each of three sites can give, then three 5-year experiments: control, pump (water removed at the sites) and inject (water added). | The author's coupled SHAKTI-ISSM intervention runs, on the model of Ehrenfeucht and Dow (2026) |

Details that differ from the sources:

- **Helheim:** the paper also solves for ice velocity; here the velocity stays at its winter value.
- **Thwaites:** pumps and injections move 90% of the capacity of each site, which is the net inflow into the site when its water pressure is held at zero.

Helheim and Thwaites have two meshes:

| Case | `low` | `high` |
|---|---|---|
| `helheimSeasonal` | 300 m to 3 km (7,366 elements) | mesh of the paper (27,913 elements) |
| `thwaitesInterventions` | 500 m to 10 km (25,046 elements) | 79,960 elements |

## Running the tests

`run/runTestCases` runs the cases you choose, then their diagnostics. It takes the cases (`all`, one name, or several), the number of processors, a run directory, and optionally the resolution (`low` by default).

MATLAB, from the `run` folder:

```matlab
runTestCases('all',4,'~/shaktiRuns')
runTestCases('helheimSeasonal',8,'~/shaktiRuns','high')
runTestCases({'singleMoulin','thwaitesInterventions'},2,'~/shaktiRuns')
```

Python, from the `run` folder:

```
python runTestCases.py all 4 ~/shaktiRuns
python runTestCases.py helheimSeasonal 8 ~/shaktiRuns high
python runTestCases.py singleMoulin,thwaitesInterventions 2 ~/shaktiRuns
```

Nothing is written in this repository. Each case runs in `<run directory>/<case>`, where its scripts and input files are copied first. If one case fails, the others still run.

You can also run a `runme` and then a `diagnostics` script by hand in the folder of a case; the results are then written in that folder.

## Output

In `<run directory>/<case>`:

- **Results:** one model file per run (`.mat` from MATLAB, `.nc` from Python), for example `Helheim_seasonal_low` or `Thwaites_pump_high`.
- **`diagnostics/<run>.mp4`:** maps of effective pressure, water flux, gap height and Reynolds number at every saved time step.
- **`diagnostics/<run>_first.png` and `<run>_last.png`:** the same maps at the first and last saved time steps.
- **`diagnostics/<run>_timeseries.png`:** water put into or taken out of the bed, and mean effective pressure over all grounded ice and over one region of the case (downstream of the moulin, below 900 m elevation at Helheim, within 200 km of the grounding line at Thwaites).

## Requirements

- ISSM with SHAKTI, and the `ISSM_DIR` environment variable.
- For Python: ISSM built with its Python interface, a version of ISSM from mid-November 2026 or later (`write_netCDF.py` has a syntax error in versions from 15 September to 12 November 2026), and the `netCDF4`, `scipy` and `matplotlib` packages. Videos are `.mp4` if `ffmpeg` or the `imageio-ffmpeg` package is installed, and `.gif` otherwise.

## Water extraction with the Python interface

As of November 2026, the Python interface of ISSM refuses a negative `moulin_input`, so the `pump` experiment of `thwaitesInterventions` stops with "Model not consistent". The MATLAB interface accepts it.

To fix it, edit `$ISSM_DIR/src/m/classes/hydrologyshakti.py` and remove `'>=', 0, ` from the check of `moulin_input` in `checkconsistency`:

```python
# before
md = checkfield(md, 'fieldname', 'hydrology.moulin_input', '>=', 0, 'NaN', 1, 'Inf', 1, 'timeseries', 1)
# after
md = checkfield(md, 'fieldname', 'hydrology.moulin_input', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
```

No recompilation is needed.

## References

- Sommers, A., Rajaram, H., and Morlighem, M. (2018). SHAKTI: Subglacial Hydrology and Kinetic, Transient Interactions v1.0. *Geoscientific Model Development*, 11, 2955–2974. https://doi.org/10.5194/gmd-11-2955-2018
- Sommers, A., et al. (2024). *Geophysical Research Letters*. https://doi.org/10.1029/2024GL109168
