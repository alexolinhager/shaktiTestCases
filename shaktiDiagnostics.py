"""Video, stills and time series of a SHAKTI run (Python version of shaktiDiagnostics.m)."""
import os
import numpy as np
from matplotlib import animation
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, LogNorm, Normalize
from matplotlib.figure import Figure
from matplotlib.tri import Triangulation

NU = 1.787e-6  # kinematic viscosity of water (m^2/s), as in SHAKTI
SEQUENTIAL = np.array([[234, 242, 251], [143, 190, 240], [42, 120, 214], [20, 70, 140], [8, 30, 70]]) / 255
NEGATIVE = np.array([[226, 226, 222], [240, 160, 120], [235, 104, 52], [150, 50, 15], [80, 20, 0]]) / 255
POSITIVE = np.array([[226, 226, 222], [143, 190, 240], [42, 120, 214], [20, 70, 140], [8, 30, 70]]) / 255
COLORS = np.array([[42, 120, 214], [235, 104, 52], [27, 175, 122], [237, 161, 0]]) / 255


def shaktiDiagnostics(md, name=None, folder='diagnostics', region=None, regionname='Region', sites=None, series=None):
    """SHAKTIDIAGNOSTICS - video, stills and time series of a SHAKTI run

    Files written (in ./diagnostics by default):
       <name>.mp4             effective pressure, water flux, gap height and
                              Reynolds number at every saved time step (a .gif
                              if ffmpeg is not available)
       <name>_first.png       the same maps at the first saved time step
       <name>_last.png        the same maps at the last saved time step
       <name>_timeseries.png  water put into (or taken out of) the bed, and
                              mean effective pressure

    Only grounded ice is shown and averaged. The Reynolds number is the saved
    one if 'HydrologyReynolds' was among the requested outputs, otherwise it
    is the water flux divided by the kinematic viscosity of water (which is
    how SHAKTI computes it).

    Usage:
       out = shaktiDiagnostics(md, name='singleMoulin')

    Options:
       name        prefix of the files (default: md.miscellaneous.name)
       folder      where the files go (default: 'diagnostics')
       region      vertices over which the effective pressure is also
                   averaged (indices from 0, or logical mask)
       regionname  name of that region in the legend
       sites       x, y of points to mark on the maps, one per row (default:
                   vertices with a moulin input, if there are 100 or fewer)
       series      extra curves for the water plot: list of dictionaries with
                   keys name, time (yr) and value (m^3/s)
    """
    if name is None:
        name = md.miscellaneous.name
    os.makedirs(folder, exist_ok=True)

    nv = md.mesh.numberofvertices
    ne = md.mesh.numberofelements
    index = np.asarray(md.mesh.elements, dtype=int) - 1
    x = np.asarray(md.mesh.x, dtype=float).flatten()
    y = np.asarray(md.mesh.y, dtype=float).flatten()
    yts = md.constants.yts

    # Saved time steps (a run that stopped early can end with an empty one)
    T = [s for s in md.results.TransientSolution.steps
         if np.size(getattr(s, 'time', [])) == 1 and np.size(getattr(s, 'EffectivePressure', [])) == nv
         and np.size(getattr(s, 'HydrologyBasalFlux', [])) == ne and np.size(getattr(s, 'HydrologyGapHeight', [])) == ne]
    if len(T) == 0:
        raise RuntimeError('shaktiDiagnostics: no hydrology results in md.results.TransientSolution')
    time = np.array([float(s.time) for s in T])
    nt = time.size
    N = np.column_stack([np.asarray(s.EffectivePressure, dtype=float).flatten() for s in T]) / 1e6  # MPa
    q = np.column_stack([np.asarray(s.HydrologyBasalFlux, dtype=float).flatten() for s in T])       # m^2/s
    gap = np.column_stack([np.asarray(s.HydrologyGapHeight, dtype=float).flatten() for s in T])     # m
    if all(np.size(getattr(s, 'HydrologyReynolds', [])) == ne for s in T):
        reynolds = np.column_stack([np.asarray(s.HydrologyReynolds, dtype=float).flatten() for s in T])
    else:
        reynolds = q / NU

    # Grounded ice: the elements shown, and the area that each vertex stands for
    ice = (np.asarray(md.mask.ice_levelset).flatten() < 0) & (np.asarray(md.mask.ocean_levelset).flatten() > 0)
    active = np.all(ice[index], axis=1)
    if not np.any(active):
        raise RuntimeError('shaktiDiagnostics: no grounded ice')
    x1, x2, x3 = x[index[:, 0]], x[index[:, 1]], x[index[:, 2]]
    y1, y2, y3 = y[index[:, 0]], y[index[:, 1]], y[index[:, 2]]
    areas = 0.5 * np.abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    weights = np.bincount(index[active, :].flatten(), weights=np.repeat(areas[active] / 3, 3), minlength=nv)

    # Time in days for a run of up to 2 years, in years otherwise
    if time[-1] - time[0] <= 2:
        timescale, timename, timelabel = yts / 86400, 'day', 'Time (days)'
    else:
        timescale, timename, timelabel = 1, 'year', 'Time (yr)'
    t = time * timescale

    # 1. Water put into the bed and taken out of it (m^3/s), at 500 times from
    # the start of the run (the forcing can vary faster than results are saved)
    tin = np.linspace(min(md.timestepping.start_time, time[0]), time[-1], 500)
    moulin = attimes(md.hydrology.moulin_input, nv, tin)
    englacial = attimes(md.hydrology.englacial_input, nv, tin)  # m/yr
    curves = [{'name': 'Moulin input', 'time': tin, 'value': np.sum(np.maximum(moulin, 0), axis=0)},
              {'name': 'Englacial input', 'time': tin, 'value': weights @ np.maximum(englacial, 0) / yts},
              {'name': 'Pump output', 'time': tin, 'value': -np.sum(np.minimum(moulin, 0), axis=0)}]
    for extra in (series or []):
        curves.append({'name': extra['name'], 'time': np.asarray(extra['time'], dtype=float).flatten(),
                       'value': np.asarray(extra['value'], dtype=float).flatten()})
    curves = [c for c in curves if np.any(c['value'] != 0)]

    # 2. Mean effective pressure (MPa), weighted by area
    Nmean = weights @ N / np.sum(weights)
    Nregion = None
    if region is not None:
        inregion = np.zeros(nv, dtype=bool)
        inregion[region] = True
        wr = weights * inregion
        if np.sum(wr) == 0:
            raise RuntimeError("shaktiDiagnostics: 'region' holds no grounded ice")
        Nregion = wr @ N / np.sum(wr)

    # Points to mark on the maps
    if sites is None:
        pos = np.where(np.any(moulin != 0, axis=1))[0]
        sites = np.column_stack((x[pos], y[pos])) if pos.size <= 100 else np.zeros((0, 2))
    sites = np.atleast_2d(np.asarray(sites, dtype=float))

    # 3. Maps: same color limits for all time steps (0.1 to 99.9th percentile)
    fields = [{'name': 'Effective pressure (MPa)', 'data': N, 'shown': N[ice, :], 'log': False},
              {'name': 'Water flux (m$^2$ s$^{-1}$)', 'data': q[active, :], 'shown': q[active, :], 'log': True},
              {'name': 'Gap height (m)', 'data': gap[active, :], 'shown': gap[active, :], 'log': True},
              {'name': 'Reynolds number', 'data': reynolds[active, :], 'shown': reynolds[active, :], 'log': True}]
    if np.max(x[ice]) - np.min(x[ice]) > 5e3:
        scale, unit = 1e-3, 'km'
    else:
        scale, unit = 1, 'm'
    fig = Figure(figsize=(12.8, 9.6), dpi=100, facecolor='w', constrained_layout=True)
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 2).flatten()
    shown = Triangulation(x * scale, y * scale, index[active, :])
    margin = 0.02 * max(np.ptp(x[ice]), np.ptp(y[ice])) * scale
    handles = []
    for ax, field in zip(axes, fields):
        if np.any(~active):
            ax.tripcolor(Triangulation(x * scale, y * scale, index[~active, :]), facecolors=np.zeros(np.sum(~active)),
                         cmap=ListedColormap([[0.9, 0.9, 0.9]]))
        values = field['shown'].flatten()
        if field['log']:
            values = values[values > 0]
            if values.size == 0:
                values = np.array([1.])
        limits = np.percentile(values, [0.1, 99.9])
        if field['log']:
            if limits[1] < 100 * limits[0]:  # whole decades, so that the colorbar has labels
                limits = 10.**np.array([np.floor(np.log10(limits[0])), np.ceil(np.log10(limits[1]))])
                if limits[1] <= limits[0]:
                    limits[1] = 10 * limits[0]
            norm = LogNorm(limits[0], limits[1], clip=True)
        else:
            if limits[1] <= limits[0]:
                limits = limits[0] + np.array([-0.5, 0.5]) * max(abs(limits[0]), 1e-6)
            norm = Normalize(limits[0], limits[1])
        if not field['log'] and limits[0] < 0 < limits[1]:
            cmap = diverging(256, limits)  # water pressure above overburden stands out
        else:
            cmap = LinearSegmentedColormap.from_list('sequential', SEQUENTIAL, 256)
        if field['log']:   # on elements
            handle = ax.tripcolor(shown, facecolors=np.maximum(field['data'][:, 0], limits[0]), cmap=cmap, norm=norm)
        else:              # on vertices
            handle = ax.tripcolor(shown, field['data'][:, 0], shading='gouraud', cmap=cmap, norm=norm)
        handles.append((handle, limits))
        if sites.size > 0:
            ax.plot(sites[:, 0] * scale, sites[:, 1] * scale, linestyle='none', marker='o', markersize=7,
                    markerfacecolor='none', markeredgecolor='k', markeredgewidth=1)
        ax.set_aspect('equal')
        ax.set_xlim(np.min(x[ice]) * scale - margin, np.max(x[ice]) * scale + margin)
        ax.set_ylim(np.min(y[ice]) * scale - margin, np.max(y[ice]) * scale + margin)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(direction='out', labelsize=11)
        ax.set_xlabel('x (' + unit + ')', fontsize=11)
        ax.set_ylabel('y (' + unit + ')', fontsize=11)
        ax.set_title(field['name'], fontsize=13)
        fig.colorbar(handle, ax=ax)

    # Video, with stills of the first and last time steps
    ffmpeg = animation.writers.is_available('ffmpeg')
    if not ffmpeg:
        try:  # the imageio-ffmpeg package comes with its own ffmpeg
            import imageio_ffmpeg
            import matplotlib
            matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
            ffmpeg = True
        except ImportError:
            pass
    fps = int(min(max(round(nt / 20), 2), 30))
    if ffmpeg:
        writer, videofile = animation.FFMpegWriter(fps=fps, bitrate=4000), os.path.join(folder, name + '.mp4')
    else:
        writer, videofile = animation.PillowWriter(fps=fps), os.path.join(folder, name + '.gif')
    title = fig.suptitle('', fontsize=15)
    with writer.saving(fig, videofile, dpi=100):
        for k in range(nt):
            for (handle, limits), field in zip(handles, fields):
                handle.set_array(np.maximum(field['data'][:, k], limits[0]) if field['log'] else field['data'][:, k])
            title.set_text('%s: %s %.2f' % (name, timename, t[k]))
            writer.grab_frame()
            if k == 0:
                fig.savefig(os.path.join(folder, name + '_first.png'), dpi=150)
            if k == nt - 1:
                fig.savefig(os.path.join(folder, name + '_last.png'), dpi=150)

    # 4. Time series: water on top (if there is any), mean effective pressure below
    nrows = 1 + (len(curves) > 0)
    fig = Figure(figsize=(9, 3 * nrows), dpi=100, facecolor='w', constrained_layout=True)
    FigureCanvasAgg(fig)
    axes = np.atleast_1d(fig.subplots(nrows, 1))
    if len(curves) > 0:
        ax = axes[0]
        for i, c in enumerate(curves):
            ax.plot(c['time'] * timescale, c['value'], linewidth=2, color=COLORS[i % 4], label=c['name'])
        ax.set_ylabel('Water (m$^3$ s$^{-1}$)', fontsize=11)
        if len(curves) > 1:
            ax.legend(frameon=False)
        else:
            ax.set_title(curves[0]['name'], fontsize=13)
    ax = axes[-1]
    ax.plot(t, Nmean, linewidth=2, color=COLORS[0], label='All grounded ice')
    if Nregion is not None:
        ax.plot(t, Nregion, linewidth=2, color=COLORS[1], label=regionname)
        ax.legend(frameon=False)
    ax.set_ylabel('Mean effective pressure (MPa)', fontsize=11)
    ax.set_xlabel(timelabel, fontsize=11)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(direction='out', labelsize=11)
        ax.grid(True, color=[0.5, 0.5, 0.5], alpha=0.25)
        if time[-1] > tin[0]:
            ax.set_xlim(tin[0] * timescale, time[-1] * timescale)
    fig.suptitle(name, fontsize=15)
    fig.savefig(os.path.join(folder, name + '_timeseries.png'), dpi=150)

    print('   diagnostics of ' + name + ' (' + str(nt) + ' time steps) written in ' + folder)
    return {'time': time, 'series': curves, 'effective_pressure': Nmean, 'effective_pressure_region': Nregion}


def attimes(field, n, time):
    # Field on vertices at the given times: one column per time. A transient
    # field holds its time stamps in its last row; ISSM interpolates linearly
    # between them, and holds the first and last values outside
    field = np.asarray(field, dtype=float)
    if field.ndim == 2 and field.shape[0] == n + 1 and field.shape[1] > 1:
        stamps = field[-1, :]
        tq = np.clip(time, stamps[0], stamps[-1])
        k = np.clip(np.searchsorted(stamps, tq, side='right') - 1, 0, stamps.size - 2)
        f = (tq - stamps[k]) / (stamps[k + 1] - stamps[k])
        return field[:n, k] * (1 - f) + field[:n, k + 1] * f
    return np.tile(field.reshape(field.shape[0], -1)[:n, :1], (1, np.size(time)))


def diverging(n, limits):
    # Two hues around a neutral gray at zero, with the same color for the same
    # distance from zero on both sides
    values = np.linspace(limits[0], limits[1], n) / np.max(np.abs(limits))
    colors = np.zeros((n, 3))
    for j in range(3):
        colors[:, j] = np.where(values < 0,
                                np.interp(-values, np.linspace(0, 1, NEGATIVE.shape[0]), NEGATIVE[:, j]),
                                np.interp(values, np.linspace(0, 1, POSITIVE.shape[0]), POSITIVE[:, j]))
    return ListedColormap(colors)
