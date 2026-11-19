# Standalone SHAKTI simulations of Thwaites Glacier with water pumped out of, or
# injected into, the subglacial system. Slimmed down from the coupled
# SHAKTI-ISSM runs of the InterventionModeling repository
# (run/runme_coupled_thwaites.m).
#
# The ice does not evolve here: geometry, grounding line and sliding velocity
# stay as in the initial condition file.
#
# Step 1: Spin-up of the subglacial system (5 years)
# Step 2: Capacity of each intervention site (5 years): the most water that a
#         pump can remove there
# Step 3: Control, pump and inject experiments (1 year each)
#
# Input files:
#    Thwaites_SHAKTI_init.nc    mesh, mask, geometry, ice velocity, friction and
#                               geothermal flux (see md.miscellaneous.notes)
#    Thwaites_SHAKTI_init_500mto10km.nc   the same on a coarser mesh
#    TIS_island.exp             outline of the island that pins the ice shelf
#
# Output files end with the resolution: ThwaitesSpinup_<resolution>.nc,
# ThwaitesCapacity_<resolution>.nc, Thwaites_<experiment>_<resolution>.nc
#
# Python version of runme.m. Vertices and elements are numbered from 0 in the
# footprints and in the capacity file (md.mesh.elements still counts from 1)
import os
import pickle
import numpy as np
from bcgslbjacobioptions import bcgslbjacobioptions
from ContourToNodes import ContourToNodes
from generic import generic
from hydrologyshakti import hydrologyshakti
from read_netCDF import read_netCDF
from results import results
from solve import solve
from write_netCDF import write_netCDF

steps = [1, 2, 3]

# Experiments to run one after the other, each one starts from the end of the
# spin-up and is saved in Thwaites_<experiment>_<resolution>.nc:
#    'control'  no intervention
#    'pump'     water removed at the sites below
#    'inject'   water added at the sites below
experiments = ['control', 'pump', 'inject']

# Intervention sites: x, y (m), each one covers all the vertices within
# pumpRadius (m)
pumpXY = np.array([[-1429020, -441873],
                   [-1438900, -424922],
                   [-1419150, -463010]])
pumpRadius = 5e3

# Number of processors and mesh resolution: change them according to your
# machine, or set nproc and resolution before running this script (as
# run/runTestCases.py does)
#    'low'   500 m to 10 km (12,690 vertices)
#    'high'  mesh of the InterventionModeling runs (40,167 vertices)
if 'nproc' not in globals():
    nproc = 2
if 'resolution' not in globals():
    resolution = 'low'

# Step 2: the head is brought down to the bed over each site (or to a water
# pressure of pumpPinFrac x overburden if pumpPinFrac > 0) over pumpRampYears,
# then held until hydroPumpTestYears. The capacity of a site is the net inflow
# into it, averaged over the last year
pumpPinFrac = 0
pumpRampYears = 1
hydroPumpTestYears = 5

# Step 3: each site removes (pump) or receives (inject) pumpFromCapacity x its
# capacity, shared between its vertices like the inflow of step 2 and ramped up
# over rampYears. A capacity above pumpCapacityMax (m^3/s) means that step 2
# had not converged
pumpFromCapacity = 0.9
pumpCapacityMax = 20
rampYears = 1


# ======================= functions (from InterventionModeling, SHAKTI branches only) =======================
def patchFootprint(md, xy, radius, pinned):
    # Vertices within radius of a site (or the nearest one): all of them (cand) and the unpinned ones (v).
    d = np.hypot(md.mesh.x - xy[0], md.mesh.y - xy[1])
    if radius > 0:
        cand = np.where(d <= radius)[0]
    else:
        cand = np.array([np.argmin(d)])
    v = cand[~pinned[cand]]
    return v, cand


def pumpPins(md, vtx, t0, tRamp, tEnd, pinFrac=0):
    # Dirichlet pins that make a pump out of each footprint: the head at its vertices is ramped linearly from the
    # current value to the target over tRamp, then held to tEnd; existing pins are kept.
    # Target: pinFrac = 0 -> the bed (p_w = 0); 0 < pinFrac < 1 -> p_w = pinFrac x overburden, i.e. head = base +
    # pinFrac rho_i/rho_fw H.
    # The pump then removes exactly the water that reaches the footprint, so it can never overdraw. ISSM relaxes a
    # transient constraint wherever its interpolated value is NaN, so the footprint must be finite from t0 on: it
    # starts at the current state.
    assert 0 <= pinFrac < 1, 'pumpPins: pinFrac must be in [0, 1)'
    nv = md.mesh.numberofvertices
    v = np.unique(np.concatenate(vtx))
    target = md.geometry.base + pinFrac * md.materials.rho_ice / md.materials.rho_freshwater * md.geometry.thickness
    spc0 = np.array(md.hydrology.spchead, dtype=float)
    cur = np.array(md.hydrology.head, dtype=float)
    assert spc0.shape == (nv,), 'pumpPins: the existing pins must be a static field on the vertices'
    assert np.all(np.isnan(spc0[v])), 'pumpPins: a footprint vertex is already pinned'
    assert cur.shape == (nv,), 'pumpPins: no current hydrology state to ramp from'
    S = np.tile(spc0[:, None], (1, 3))
    S[v, 0] = cur[v]
    S[v, 1] = target[v]
    S[v, 2] = target[v]
    times = np.array([t0, t0 + tRamp, tEnd])
    if tRamp <= 0:
        S = S[:, [1, 2]]
        times = np.array([t0, tEnd])
    md.hydrology.spchead = np.vstack((S, times))
    md.hydrology.moulin_input = np.zeros(nv)
    md.hydrology.englacial_input = np.zeros(nv)
    return md


def elementFlux(md, S):
    # Sheet flux vector of one output step S, rebuilt per element from the saved fields the way the solver forms it:
    #   SHAKTI  q = HydrologyBasalFlux * (-grad h / |grad h|)         (BasalFlux = K |grad h| per element)
    # Elements that are all floating carry no flux. Also returns the P1 gradient coefficients (df/dx =
    # sum(bx * f[el], axis=1), df/dy = sum(by * f[el], axis=1)) and the element areas.
    el = md.mesh.elements.astype(int) - 1
    x = md.mesh.x
    y = md.mesh.y
    x1, x2, x3 = x[el[:, 0]], x[el[:, 1]], x[el[:, 2]]
    y1, y2, y3 = y[el[:, 0]], y[el[:, 1]], y[el[:, 2]]
    twoA = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    bx = np.column_stack((y2 - y3, y3 - y1, y1 - y2)) / twoA[:, None]
    by = np.column_stack((x3 - x2, x1 - x3, x2 - x1)) / twoA[:, None]
    f = np.asarray(S.HydrologyHead, dtype=float).flatten()
    gx = np.sum(bx * f[el], axis=1)
    gy = np.sum(by * f[el], axis=1)
    gn = np.maximum(np.hypot(gx, gy), np.finfo(float).eps)
    coeff = np.asarray(S.HydrologyBasalFlux, dtype=float).flatten() / gn
    coeff[np.all(md.mask.ocean_levelset[el] <= 0, axis=1)] = 0   # the ice does not evolve: the mask of the model holds
    return -coeff * gx, -coeff * gy, bx, by, np.abs(twoA) / 2


def nodalInflow(md, S, v):
    # Water arriving at each vertex of v (m^3/s) from the sheet flux of one output step S: the discrete divergence
    # sum_e A_e (q_e . grad(phi_i)).
    # Used for the share weights of a prescribed discharge, so only the distribution matters.
    el = md.mesh.elements.astype(int) - 1
    qx, qy, bx, by, A = elementFlux(md, S)
    contrib = A[:, None] * (qx[:, None] * bx + qy[:, None] * by)   # ne x 3: flux into each vertex of each element
    qn = np.bincount(el.flatten(), weights=contrib.flatten(), minlength=md.mesh.numberofvertices)
    return qn[v]


def patchWaterFlux(md, T, v):
    # Water flux (m^3/s) leaving and entering the footprint of vertices v at every saved step of T: the sheet flux
    # integrated as q.n over the perimeter of the footprint elements (elements with all vertices in v).
    el = md.mesh.elements.astype(int) - 1
    x = md.mesh.x
    y = md.mesh.y
    ne = el.shape[0]
    inV = np.zeros(md.mesh.numberofvertices, dtype=bool)
    inV[v] = True
    elIn = np.all(inV[el], axis=1)
    if not np.any(elIn):
        elIn = np.any(inV[el], axis=1)                               # tiny footprint: elements touching it

    # unique edges and their one or two adjacent elements (-1: none, the edge is on the domain boundary)
    pairs = np.sort(np.vstack((el[:, [0, 1]], el[:, [1, 2]], el[:, [2, 0]])), axis=1)
    edgeV, idx = np.unique(pairs, axis=0, return_inverse=True)
    idx = idx.flatten()
    eid = np.tile(np.arange(ne), 3)
    e1 = np.full(edgeV.shape[0], ne)
    np.minimum.at(e1, idx, eid)
    e2 = np.full(edgeV.shape[0], -1)
    np.maximum.at(e2, idx, eid)
    e2[np.bincount(idx, minlength=edgeV.shape[0]) == 1] = -1
    side1 = elIn[e1]
    side2 = (e2 >= 0) & elIn[np.maximum(e2, 0)]
    perim = np.where(side1 ^ side2)[0]                               # exactly one adjacent element inside
    inner = np.where(side1[perim], e1[perim], e2[perim])
    outer = np.where(side1[perim], e2[perim], e1[perim])
    outer = np.where(outer < 0, inner, outer)                        # boundary edge: use the inside element

    # outward unit normal and length of each perimeter edge
    ex = x[edgeV[perim, :]]
    ey = y[edgeV[perim, :]]
    L = np.hypot(ex[:, 1] - ex[:, 0], ey[:, 1] - ey[:, 0])
    nx = (ey[:, 1] - ey[:, 0]) / L
    ny = -(ex[:, 1] - ex[:, 0]) / L
    cx = np.mean(x[el[inner, :]], axis=1)
    cy = np.mean(y[el[inner, :]], axis=1)
    flip = nx * (np.mean(ex, axis=1) - cx) + ny * (np.mean(ey, axis=1) - cy) < 0
    nx[flip] = -nx[flip]
    ny[flip] = -ny[flip]

    outflow = np.full(len(T), np.nan)
    inflow = np.full(len(T), np.nan)
    for i in range(len(T)):
        qx, qy = elementFlux(md, T[i])[0:2]
        qn = qx[outer] * nx + qy[outer] * ny
        outflow[i] = np.sum(np.maximum(0, qn) * L)
        inflow[i] = np.sum(np.maximum(0, -qn) * L)
    return outflow, inflow


def savedSteps(md):
    # Output steps of the last run that hold results (a run that stopped early can end with an empty one)
    nv = md.mesh.numberofvertices
    return [s for s in md.results.TransientSolution.steps
            if np.size(getattr(s, 'time', [])) == 1 and np.size(getattr(s, 'EffectivePressure', [])) == nv]


def restart(md):
    # Start from the end of the spin-up
    spinup = savedSteps(md)[-1]
    md.hydrology.head = np.asarray(spinup.HydrologyHead, dtype=float).flatten()
    md.hydrology.gap_height = np.asarray(spinup.HydrologyGapHeight, dtype=float).flatten()
    md.hydrology.reynolds = np.maximum(np.asarray(spinup.HydrologyBasalFlux, dtype=float).flatten() / 1.787e-6, 1)
    md.results = results()
    return md


def footprints(md):
    # every vertex within pumpRadius of each site, minus pinned ones (a Dirichlet head cannot take a flux)
    pinned = ~np.isnan(md.hydrology.spchead)
    vtx = []
    for i in range(pumpXY.shape[0]):
        vtx.append(patchFootprint(md, pumpXY[i, :], pumpRadius, pinned)[0])
        assert vtx[i].size > 0, 'pump %d: every vertex within %g m of (%.0f, %.0f) is pinned' % (i + 1, pumpRadius, pumpXY[i, 0], pumpXY[i, 1])
    return vtx


# ======================= steps =======================
if 1 in steps:
    print('	Step 1: Spin-up')
    if resolution == 'low':
        md = read_netCDF('Thwaites_SHAKTI_init_500mto10km.nc')
    elif resolution == 'high':
        md = read_netCDF('Thwaites_SHAKTI_init.nc')
    else:
        raise ValueError("resolution '" + resolution + "' not supported")
    nv = md.mesh.numberofvertices
    ne = md.mesh.numberofelements

    # Standalone SHAKTI: only solve for hydrology
    md.transient = md.transient.deactivateall()
    md.transient.ishydrology = 1
    md.hydrology = hydrologyshakti()

    rho_ice = md.materials.rho_ice
    rho_freshwater = md.materials.rho_freshwater
    rho_seawater = md.materials.rho_water

    # Define initial water head such that water pressure is 50% of ice overburden pressure
    md.hydrology.head = 0.5 * rho_ice / rho_freshwater * md.geometry.thickness + md.geometry.base

    # Initial subglacial gap height of 0.01m everywhere
    md.hydrology.gap_height = 0.01 * np.ones(ne)

    # Initial Reynolds number (start at Re=1000 everywhere)
    md.hydrology.reynolds = 1000. * np.ones(ne)

    # No opening of the gap by sliding over bed bumps (bump height of 0m)
    md.hydrology.bump_spacing = 1.0 * np.ones(ne)
    md.hydrology.bump_height = 0.0 * np.ones(ne)

    # No water input other than the melt computed by SHAKTI
    md.hydrology.englacial_input = np.zeros(nv)
    md.hydrology.moulin_input = np.zeros(nv)

    # No-flux Type 2 boundary conditions on all edges (except where the head
    # is constrained below)
    md.hydrology.neumannflux = np.zeros(ne)

    # Type 1 boundary conditions: head of the ocean under floating ice and in
    # open water, and effective pressure held at 30% of ice overburden pressure
    # on the grounded part of the island
    md.hydrology.spchead = np.nan * np.ones(nv)
    pos = np.where(md.mask.ocean_levelset <= 0)[0]
    md.hydrology.spchead[pos] = md.geometry.bed[pos] * (rho_freshwater - rho_seawater) / rho_freshwater
    island = np.asarray(ContourToNodes(md.mesh.x, md.mesh.y, 'TIS_island.exp', 2)).flatten()[:nv] > 0
    pos = np.where(island & (md.mask.ocean_levelset > 0))[0]
    md.hydrology.spchead[pos] = md.geometry.bed[pos] + (1 - 0.3) * rho_ice / rho_freshwater * md.geometry.thickness[pos]

    # Convergence criteria of the SHAKTI nonlinear solver, and iterative linear solver
    md.hydrology.maxiter = 40
    md.hydrology.restol = 0.01
    md.hydrology.reltol = 0.1
    md.hydrology.abstol = np.nan
    md.toolkits.DefaultAnalysis = bcgslbjacobioptions()
    md.settings.solver_residue_threshold = np.nan

    # Define the time stepping scheme: run for 5 years with a time step of 6 hr,
    # and save the model output every 15 days
    md.timestepping.start_time = 0
    md.timestepping.time_step = 6 * 3600 / md.constants.yts  # Time step (in years)
    md.timestepping.final_time = 5
    md.settings.output_frequency = 2 * 30

    # Specify that you want to run the model on your current computer
    md.cluster = generic('np', nproc)

    md.miscellaneous.name = 'ThwaitesSpinup_' + resolution
    md.verbose.solution = True
    md = solve(md, 'Transient')

    write_netCDF(md, 'ThwaitesSpinup_' + resolution + '.nc')

if 2 in steps:
    print('	Step 2: Capacity of each site')
    md = read_netCDF('ThwaitesSpinup_' + resolution + '.nc')
    md = restart(md)
    nv = md.mesh.numberofvertices

    # A. footprints
    vtx = footprints(md)
    nP = len(vtx)
    allv = np.concatenate(vtx)
    assert np.unique(allv).size == allv.size, 'step 2: footprints overlap; reduce pumpRadius or move the sites'

    # B. hydrology on the state and settings of the experiments, footprints pinned to the bed
    md.timestepping.final_time = hydroPumpTestYears
    md.cluster = generic('np', nproc)
    md = pumpPins(md, vtx, 0, pumpRampYears, hydroPumpTestYears, pumpPinFrac)
    md.miscellaneous.name = 'ThwaitesCapacity_' + resolution
    md = solve(md, 'Transient')
    write_netCDF(md, 'ThwaitesCapacity_' + resolution + '.nc')
    # a solver failure returns whatever was saved before it; a capacity from that is meaningless, so stop here
    T = savedSteps(md)
    tReached = max([float(s.time) for s in T] + [0])
    assert tReached >= hydroPumpTestYears - 2 * md.timestepping.time_step - 1e-9, \
        'step 2: the capacity run stopped at %.2f of %g yr (solver failure, see the log); no capacity written. Lower the drawdown (pumpPinFrac > 0), the number of sites, or the time step' % (tReached, hydroPumpTestYears)

    # C. capacity = net inflow into each pinned footprint, averaged over the last year; its range over that year shows
    #    whether the solution settled (SHAKTI can flip between two states next to a dry-bed pin); the per-vertex share
    #    of the inflow is the weight step 3 uses to split a prescribed discharge
    time = np.array([float(s.time) for s in T])
    late = time >= hydroPumpTestYears - 1 - 1e-9
    if not np.any(late):
        late[-1] = True
    cap = {'source': md.miscellaneous.name, 'xy': pumpXY, 'vertices': vtx, 'radius': pumpRadius, 'pinFrac': pumpPinFrac,
           'time': time, 'Qin': np.zeros((nP, len(T))), 'capacity': np.zeros(nP), 'capacityMin': np.zeros(nP),
           'capacityMax': np.zeros(nP), 'weights': [None] * nP}
    for i in range(nP):
        o, f = patchWaterFlux(md, T, vtx[i])
        cap['Qin'][i, :] = f - o
        cap['capacity'][i] = np.mean(cap['Qin'][i, late])
        cap['capacityMin'][i] = np.min(cap['Qin'][i, late])
        cap['capacityMax'][i] = np.max(cap['Qin'][i, late])
        w = np.mean([nodalInflow(md, T[k], vtx[i]) for k in np.where(late)[0]], axis=0)   # steps x vertices -> mean per vertex
        w = np.maximum(w, 0)
        if np.sum(w) > 0:
            w = w / np.sum(w)
        else:
            w = np.ones(vtx[i].size) / vtx[i].size
        cap['weights'][i] = w
    with open('ThwaitesCapacity_cap_' + resolution + '.pkl', 'wb') as fid:
        pickle.dump(cap, fid)

    # D. report
    print('\nstep 2: pump capacity, %g yr (head held at pumpPinFrac = %g over each footprint after a %g yr ramp)' % (hydroPumpTestYears, pumpPinFrac, pumpRampYears))
    print('%-6s %-6s | %-9s %-9s %-9s' % ('pump', 'nvert', 'capacity', 'min', 'max'))
    for i in range(nP):
        print('%-6d %-6d | %-9.3g %-9.3g %-9.3g' % (i + 1, vtx[i].size, cap['capacity'][i], cap['capacityMin'][i], cap['capacityMax'][i]))
    print('capacity in m^3/s (net inflow into the pinned footprint); a wide min-max range means the solution has not settled')

if 3 in steps:
    for experiment in experiments:
        print('	Step 3: Experiment (' + experiment + ')')
        md = read_netCDF('ThwaitesSpinup_' + resolution + '.nc')
        md = restart(md)
        nv = md.mesh.numberofvertices

        # Run for 5 years
        md.timestepping.final_time = 5
        t0 = md.timestepping.start_time
        finalYear = md.timestepping.final_time

        if experiment == 'control':
            factor = 0
        elif experiment == 'pump':
            factor = -1
        elif experiment == 'inject':
            factor = 1
        else:
            raise ValueError("experiment '" + experiment + "' not supported")

        if factor != 0:
            assert rampYears < finalYear - t0, 'rampYears must be shorter than the run (%g yr)' % (finalYear - t0)

            # A. footprints
            vtx = footprints(md)
            nP = len(vtx)

            # B. discharges and per-vertex shares: a fraction of the step 2 capacity split by the capacity weights
            capFile = 'ThwaitesCapacity_cap_' + resolution + '.pkl'
            assert os.path.isfile(capFile), 'step 3: needs the step 2 capacity file %s' % capFile
            with open(capFile, 'rb') as fid:
                cap = pickle.load(fid)
            assert np.array_equal(cap['xy'], pumpXY) and cap['radius'] == pumpRadius and len(cap['vertices']) == nP \
                and all(np.array_equal(a, b) for a, b in zip(cap['vertices'], vtx)), \
                'step 3: %s was made for other sites or footprints; rerun step 2' % capFile
            assert cap['pinFrac'] == pumpPinFrac, 'step 3: %s was made with pumpPinFrac = %g, not %g' % (capFile, cap['pinFrac'], pumpPinFrac)
            # sanity: a diverged capacity run gives huge or non-finite values; a slightly negative one is a footprint
            # that exported water (nothing to pump there)
            assert np.all(np.isfinite(cap['capacity'])), 'step 3: non-finite capacity in %s; rerun step 2' % capFile
            assert np.max(cap['capacity']) <= pumpCapacityMax, \
                'step 3: capacity up to %.3g m^3/s in %s exceeds pumpCapacityMax = %g: the step 2 run had not converged; rerun it (pumpPinFrac > 0 near the grounding line)' % (np.max(cap['capacity']), capFile, pumpCapacityMax)
            if np.any(cap['capacity'] < 0):
                print('step 3: %d sites with a negative capacity (footprint exported water) set to 0' % np.sum(cap['capacity'] < 0))
                cap['capacity'] = np.maximum(cap['capacity'], 0)
            Qp = pumpFromCapacity * cap['capacity']
            wts = cap['weights']
            for i in range(nP):
                print('%s %d: %.3g m^3/s = %g x capacity %.3g (min %.3g over the last year of step 2) over %d vertices (max share %.3g)'
                      % (experiment, i + 1, Qp[i], pumpFromCapacity, cap['capacity'][i], cap['capacityMin'][i], vtx[i].size, Qp[i] * np.max(wts[i])))

            # C. per-vertex input at full strength (m^3/s per vertex, additive where sites overlap) and the series:
            #    0 -> full over rampYears, then held (ISSM interpolates linearly between stamps)
            Mv = np.zeros(nv)
            for i in range(nP):
                np.add.at(Mv, vtx[i], Qp[i] * wts[i])
            md.hydrology.moulin_input = np.vstack((factor * np.column_stack((np.zeros(nv), Mv, Mv)),
                                                   [t0, t0 + rampYears, finalYear]))   # last row = time stamps (yr)
            print('step 3 (%s): %d sites, %.3g m^3/s in total, ramped over %g yr, held to %g yr' % (experiment, nP, np.sum(Qp), rampYears, finalYear))

        # Define the cluster again, in case the name of your computer has changed
        # since the spin-up
        md.cluster = generic('np', nproc)

        md.miscellaneous.name = 'Thwaites_' + experiment + '_' + resolution
        md = solve(md, 'Transient')

        write_netCDF(md, 'Thwaites_' + experiment + '_' + resolution + '.nc')
