% Standalone SHAKTI simulations of Thwaites Glacier with water pumped out of, or
% injected into, the subglacial system. Slimmed down from the coupled
% SHAKTI-ISSM runs of the InterventionModeling repository
% (run/runme_coupled_thwaites.m).
%
% The ice does not evolve here: geometry, grounding line and sliding velocity
% stay as in the initial condition file.
%
% Step 1: Spin-up of the subglacial system (5 years)
% Step 2: Capacity of each intervention site (5 years): the most water that a
%         pump can remove there
% Step 3: Control, pump and inject experiments (1 year each)
%
% Input files:
%    Thwaites_SHAKTI_init.mat   mesh, mask, geometry, ice velocity, friction and
%                               geothermal flux (see md.miscellaneous.notes)
%    Thwaites_SHAKTI_init_500mto10km.mat   the same on a coarser mesh
%    TIS_island.exp             outline of the island that pins the ice shelf
%
% Output files end with the resolution: ThwaitesSpinup_<resolution>.mat,
% ThwaitesCapacity_<resolution>.mat, Thwaites_<experiment>_<resolution>.mat

steps=[1:3];

% Experiments to run one after the other, each one starts from the end of the
% spin-up and is saved in Thwaites_<experiment>_<resolution>.mat:
%    'control'  no intervention
%    'pump'     water removed at the sites below
%    'inject'   water added at the sites below
experiments={'control','pump','inject'};

% Intervention sites: x, y (m), each one covers all the vertices within
% pumpRadius (m)
pumpXY=[-1429020 -441873
        -1438900 -424922
        -1419150 -463010];
pumpRadius=5e3;

% Number of processors and mesh resolution: change them according to your
% machine, or set np and resolution before running this script (as
% run/runTestCases.m does)
%    'low'   500 m to 10 km (12,690 vertices)
%    'high'  mesh of the InterventionModeling runs (40,167 vertices)
if ~exist('np','var'), np=2; end
if ~exist('resolution','var'), resolution='low'; end

% Step 2: the head is brought down to the bed over each site (or to a water
% pressure of pumpPinFrac x overburden if pumpPinFrac > 0) over pumpRampYears,
% then held until hydroPumpTestYears. The capacity of a site is the net inflow
% into it, averaged over the last year
pumpPinFrac=0;
pumpRampYears=1;
hydroPumpTestYears=5;

% Step 3: each site removes (pump) or receives (inject) pumpFromCapacity x its
% capacity, shared between its vertices like the inflow of step 2 and ramped up
% over rampYears. A capacity above pumpCapacityMax (m^3/s) means that step 2
% had not converged
pumpFromCapacity=0.9;
pumpCapacityMax=20;
rampYears=1;

if any(steps==1)
	disp('	Step 1: Spin-up');
	switch resolution
		case 'low',  md=loadmodel('Thwaites_SHAKTI_init_500mto10km');
		case 'high', md=loadmodel('Thwaites_SHAKTI_init');
		otherwise, error(['resolution ''' resolution ''' not supported']);
	end

	% Standalone SHAKTI: only solve for hydrology
	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;
	md.hydrology=hydrologyshakti();

	rho_ice=md.materials.rho_ice;
	rho_freshwater=md.materials.rho_freshwater;
	rho_seawater=md.materials.rho_water;

	% Define initial water head such that water pressure is 50% of ice overburden pressure
	md.hydrology.head=0.5*rho_ice/rho_freshwater*md.geometry.thickness+md.geometry.base;

	% Initial subglacial gap height of 0.01m everywhere
	md.hydrology.gap_height=0.01*ones(md.mesh.numberofelements,1);

	% Initial Reynolds number (start at Re=1000 everywhere)
	md.hydrology.reynolds=1000*ones(md.mesh.numberofelements,1);

	% No opening of the gap by sliding over bed bumps (bump height of 0m)
	md.hydrology.bump_spacing=1.0*ones(md.mesh.numberofelements,1);
	md.hydrology.bump_height=0.0*ones(md.mesh.numberofelements,1);

	% No water input other than the melt computed by SHAKTI
	md.hydrology.englacial_input=zeros(md.mesh.numberofvertices,1);
	md.hydrology.moulin_input=zeros(md.mesh.numberofvertices,1);

	% No-flux Type 2 boundary conditions on all edges (except where the head
	% is constrained below)
	md.hydrology.neumannflux=zeros(md.mesh.numberofelements,1);

	% Type 1 boundary conditions: head of the ocean under floating ice and in
	% open water, and effective pressure held at 30% of ice overburden pressure
	% on the grounded part of the island
	md.hydrology.spchead=NaN(md.mesh.numberofvertices,1);
	pos=find(md.mask.ocean_levelset<=0);
	md.hydrology.spchead(pos)=md.geometry.bed(pos)*(rho_freshwater-rho_seawater)/rho_freshwater;
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'TIS_island.exp',2) & md.mask.ocean_levelset>0);
	md.hydrology.spchead(pos)=md.geometry.bed(pos)+(1-0.3)*rho_ice/rho_freshwater*md.geometry.thickness(pos);

	% Convergence criteria of the SHAKTI nonlinear solver, and iterative linear solver
	md.hydrology.maxiter=40;
	md.hydrology.restol=0.01;
	md.hydrology.reltol=0.1;
	md.hydrology.abstol=NaN;
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();
	md.settings.solver_residue_threshold=NaN;

	% Define the time stepping scheme: run for 10 years with a time step of 12 hr,
	% and save the model output every 30 days
	md.timestepping.start_time=0;
	md.timestepping.time_step=6*3600/md.constants.yts; % Time step (in years)
	md.timestepping.final_time=5;
	md.settings.output_frequency=2*30;

	% Specify that you want to run the model on your current computer
	md.cluster=generic('np',np);

	md.miscellaneous.name=['ThwaitesSpinup_' resolution];
	md.verbose.solution=1;
	md=solve(md,'Transient');

	save(['ThwaitesSpinup_' resolution],'md','-v7.3');
end

if any(steps==2)
	disp('	Step 2: Capacity of each site');
	md=loadmodel(['ThwaitesSpinup_' resolution]);
	nv = md.mesh.numberofvertices;

	% Start from the end of the spin-up
	spinup=md.results.TransientSolution(end);
	md.hydrology.head=spinup.HydrologyHead;
	md.hydrology.gap_height=spinup.HydrologyGapHeight;
	md.hydrology.reynolds=max(spinup.HydrologyBasalFlux./1.787e-6,1);
	md.results=struct();

	% A. footprints: every vertex within pumpRadius, minus pinned ones (a Dirichlet head cannot take a flux)
	pinned = ~isnan(md.hydrology.spchead);
	nP = size(pumpXY,1);
	vtx = cell(nP,1);
	for i = 1:nP
		vtx{i} = patchFootprint(md, pumpXY(i,:), pumpRadius, pinned);
		assert(~isempty(vtx{i}), 'pump %d: every vertex within %g m of (%.0f, %.0f) is pinned', i, pumpRadius, pumpXY(i,1), pumpXY(i,2));
	end
	allv = vertcat(vtx{:});
	assert(numel(unique(allv)) == numel(allv), 'step 2: footprints overlap; reduce pumpRadius or move the sites');

	% B. hydrology on the state and settings of the experiments, footprints pinned to the bed
	md.timestepping.final_time = hydroPumpTestYears;
	md.cluster = generic('np',np);
	md = pumpPins(md, vtx, 0, pumpRampYears, hydroPumpTestYears, pumpPinFrac);
	md.miscellaneous.name = ['ThwaitesCapacity_' resolution];
	md = solve(md,'Transient');
	save(['ThwaitesCapacity_' resolution],'md','-v7.3');
	% a solver failure returns whatever was saved before it; a capacity from that is meaningless, so stop here
	Tall = md.results.TransientSolution;
	tReached = max([Tall(arrayfun(@(st) ~isempty(st.time), Tall)).time]);
	assert(~isempty(tReached) && tReached >= hydroPumpTestYears - 2*md.timestepping.time_step - 1e-9, ...
		'step 2: the capacity run stopped at %.2f of %g yr (solver failure, see the log); no capacity written. Lower the drawdown (pumpPinFrac > 0), the number of sites, or the time step', ...
		max([0 tReached]), hydroPumpTestYears);

	% C. capacity = net inflow into each pinned footprint, averaged over the last year; its range over that year shows
	%    whether the solution settled (SHAKTI can flip between two states next to a dry-bed pin); the per-vertex share
	%    of the inflow is the weight step 3 uses to split a prescribed discharge
	T = md.results.TransientSolution;
	T = T(arrayfun(@(st) numel(st.EffectivePressure) == nv && ~isempty(st.time), T));
	T = fillIceFields(md, T);
	late = [T.time] >= hydroPumpTestYears - 1 - 1e-9;
	if ~any(late), late(end) = true; end
	cap = struct('source', md.miscellaneous.name, 'xy', pumpXY, 'vertices', {vtx}, ...
		'radius', pumpRadius, 'pinFrac', pumpPinFrac, 'time', [T.time], 'Qin', zeros(nP, numel(T)), 'capacity', zeros(nP,1), ...
		'capacityMin', zeros(nP,1), 'capacityMax', zeros(nP,1), 'weights', {cell(nP,1)});
	for i = 1:nP
		[o, f] = patchWaterFlux(md, T, vtx{i});
		cap.Qin(i,:) = f - o;
		cap.capacity(i)    = mean(cap.Qin(i, late));
		cap.capacityMin(i) = min(cap.Qin(i, late));
		cap.capacityMax(i) = max(cap.Qin(i, late));
		W = arrayfun(@(st) nodalInflow(md, st, vtx{i})', T(late), 'UniformOutput', false);
		w = mean(vertcat(W{:}), 1)';                                   % steps x vertices -> mean per vertex
		w = max(w, 0);
		if sum(w) > 0, w = w / sum(w); else, w = ones(numel(vtx{i}),1) / numel(vtx{i}); end
		cap.weights{i} = w;
	end
	save(['ThwaitesCapacity_cap_' resolution],'cap');

	% D. report
	fprintf('\nstep 2: pump capacity, %g yr (head held at pumpPinFrac = %g over each footprint after a %g yr ramp)\n', hydroPumpTestYears, pumpPinFrac, pumpRampYears);
	fprintf('%-6s %-6s | %-9s %-9s %-9s\n', 'pump', 'nvert', 'capacity', 'min', 'max');
	for i = 1:nP
		fprintf('%-6d %-6d | %-9.3g %-9.3g %-9.3g\n', i, numel(vtx{i}), cap.capacity(i), cap.capacityMin(i), cap.capacityMax(i));
	end
	fprintf('capacity in m^3/s (net inflow into the pinned footprint); a wide min-max range means the solution has not settled\n');
end

if any(steps==3)
	for k=1:numel(experiments)
		experiment=experiments{k};
		disp(['	Step 3: Experiment (' experiment ')']);
		md=loadmodel(['ThwaitesSpinup_' resolution]);
		nv = md.mesh.numberofvertices;

		% Start from the end of the spin-up
		spinup=md.results.TransientSolution(end);
		md.hydrology.head=spinup.HydrologyHead;
		md.hydrology.gap_height=spinup.HydrologyGapHeight;
		md.hydrology.reynolds=max(spinup.HydrologyBasalFlux./1.787e-6,1);
		md.results=struct();

		% Run for 1 year
		md.timestepping.final_time=5;
		t0 = md.timestepping.start_time;
		finalYear = md.timestepping.final_time;

		switch experiment
			case 'control', factor=0;
			case 'pump',    factor=-1;
			case 'inject',  factor=1;
			otherwise, error(['experiment ''' experiment ''' not supported']);
		end

		if factor~=0
			assert(rampYears < finalYear - t0, 'rampYears must be shorter than the run (%g yr)', finalYear - t0);

			% A. footprints: every vertex within pumpRadius, minus pinned ones (a Dirichlet head cannot take a flux)
			pinned = ~isnan(md.hydrology.spchead);
			nP = size(pumpXY,1);
			vtx = cell(nP,1);
			for i = 1:nP
				vtx{i} = patchFootprint(md, pumpXY(i,:), pumpRadius, pinned);
				assert(~isempty(vtx{i}), 'pump %d: every vertex within %g m of (%.0f, %.0f) is pinned', i, pumpRadius, pumpXY(i,1), pumpXY(i,2));
			end

			% B. discharges and per-vertex shares: a fraction of the step 2 capacity split by the capacity weights
			capFile = ['ThwaitesCapacity_cap_' resolution '.mat'];
			assert(isfile(capFile), 'step 3: needs the step 2 capacity file %s', capFile);
			load(capFile, 'cap');
			assert(isequal(cap.xy, pumpXY) && isequal(cap.radius, pumpRadius) && all(cellfun(@isequal, cap.vertices(:), vtx(:))), ...
				'step 3: %s was made for other sites or footprints; rerun step 2', capFile);
			assert(cap.pinFrac == pumpPinFrac, 'step 3: %s was made with pumpPinFrac = %g, not %g', capFile, cap.pinFrac, pumpPinFrac);
			% sanity: a diverged capacity run gives huge or non-finite values; a slightly negative one is a footprint
			% that exported water (nothing to pump there)
			assert(all(isfinite(cap.capacity)), 'step 3: non-finite capacity in %s; rerun step 2', capFile);
			assert(max(cap.capacity) <= pumpCapacityMax, 'step 3: capacity up to %.3g m^3/s in %s exceeds pumpCapacityMax = %g: the step 2 run had not converged; rerun it (pumpPinFrac > 0 near the grounding line)', ...
				max(cap.capacity), capFile, pumpCapacityMax);
			if any(cap.capacity < 0)
				fprintf('step 3: %d sites with a negative capacity (footprint exported water) set to 0\n', nnz(cap.capacity < 0));
				cap.capacity = max(cap.capacity, 0);
			end
			Qp  = pumpFromCapacity * cap.capacity(:);
			wts = cap.weights(:);
			for i = 1:nP
				fprintf('%s %d: %.3g m^3/s = %g x capacity %.3g (min %.3g over the last year of step 2) over %d vertices (max share %.3g)\n', ...
					experiment, i, Qp(i), pumpFromCapacity, cap.capacity(i), cap.capacityMin(i), numel(vtx{i}), Qp(i)*max(wts{i}));
			end

			% C. per-vertex input at full strength (m^3/s per vertex, additive where sites overlap) and the series:
			%    0 -> full over rampYears, then held (ISSM interpolates linearly between stamps)
			Mv = zeros(nv,1);
			for i = 1:nP
				v = vtx{i};  w = wts{i}(:);
				Mv(v) = Mv(v) + Qp(i) * w;
			end
			md.hydrology.moulin_input = [factor*[zeros(nv,1) Mv Mv]; t0 t0+rampYears finalYear];   % last row = time stamps (yr)
			fprintf('step 3 (%s): %d sites, %.3g m^3/s in total, ramped over %g yr, held to %g yr\n', experiment, nP, sum(Qp), rampYears, finalYear);
		end

		% Define the cluster again, in case the name of your computer has changed
		% since the spin-up
		md.cluster=generic('np',np);

		md.miscellaneous.name=['Thwaites_' experiment '_' resolution];
		md=solve(md,'Transient');

		save(['Thwaites_' experiment '_' resolution],'md','-v7.3');
	end
end

%% ======================= local functions (from InterventionModeling, SHAKTI branches only) =======================
function [v, cand] = patchFootprint(md, xy, radius, pinned)
	% Vertices within radius of a site (or the nearest one): all of them (cand) and the unpinned ones (v).
	d = hypot(md.mesh.x(:) - xy(1), md.mesh.y(:) - xy(2));
	[~, pos] = min(d);
	if radius > 0, cand = find(d <= radius); else, cand = pos; end
	v = cand(~pinned(cand));
end

function md = pumpPins(md, vtx, t0, tRamp, tEnd, pinFrac)
	% Dirichlet pins that make a pump out of each footprint: the head at its vertices is ramped linearly from the
	% current value to the target over tRamp, then held to tEnd; existing pins are kept.
	% Target: pinFrac = 0 -> the bed (p_w = 0); 0 < pinFrac < 1 -> p_w = pinFrac x overburden, i.e. head = base +
	% pinFrac rho_i/rho_fw H.
	% The pump then removes exactly the water that reaches the footprint, so it can never overdraw. ISSM relaxes a
	% transient constraint wherever its interpolated value is NaN, so the footprint must be finite from t0 on: it
	% starts at the current state.
	if nargin < 6, pinFrac = 0; end
	assert(pinFrac >= 0 && pinFrac < 1, 'pumpPins: pinFrac must be in [0, 1)');
	nv = md.mesh.numberofvertices;  v = unique(vertcat(vtx{:}));
	head_target = md.geometry.base(:) + pinFrac * md.materials.rho_ice/md.materials.rho_freshwater * md.geometry.thickness(:);
	spc0 = md.hydrology.spchead(:);  cur = md.hydrology.head(:);  target = head_target;
	assert(numel(spc0) == nv, 'pumpPins: the existing pins must be a static nv x 1 field');
	assert(all(isnan(spc0(v))), 'pumpPins: a footprint vertex is already pinned');
	assert(numel(cur) == nv, 'pumpPins: no current hydrology state to ramp from');
	S = repmat(spc0, 1, 3);
	S(v,1) = cur(v);  S(v,2) = target(v);  S(v,3) = target(v);
	times = [t0 t0+tRamp tEnd];
	if tRamp <= 0, S = S(:,[2 3]);  times = [t0 tEnd]; end
	S = [S; times];
	md.hydrology.spchead = S;
	md.hydrology.moulin_input = zeros(nv,1);
	md.hydrology.englacial_input = zeros(nv,1);
end

function q = nodalInflow(md, S, v)
	% Water arriving at each vertex of v (m^3/s) from the sheet flux of one output step S: the discrete divergence
	% sum_e A_e (q_e . grad(phi_i)) with the element flux rebuilt as in patchWaterFlux.
	% Used for the share weights of a prescribed discharge, so only the distribution matters.
	el = md.mesh.elements;  x = md.mesh.x(:);  y = md.mesh.y(:);  nv = numel(x);
	A = GetAreas(el, x, y);
	x1 = x(el(:,1)); x2 = x(el(:,2)); x3 = x(el(:,3));  y1 = y(el(:,1)); y2 = y(el(:,2)); y3 = y(el(:,3));
	twoA = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);
	bx = [y2-y3, y3-y1, y1-y2] ./ twoA;  by = [x3-x2, x1-x3, x2-x1] ./ twoA;
	f = double(S.HydrologyHead(:));
	gx = sum(bx.*f(el), 2);  gy = sum(by.*f(el), 2);  gn = max(hypot(gx,gy), eps);
	coeff = double(S.HydrologyBasalFlux(:)) ./ gn;
	if isfield(S, 'MaskOceanLevelset') && numel(S.MaskOceanLevelset) == nv
		coeff(all(S.MaskOceanLevelset(el) <= 0, 2)) = 0;
	end
	qx = -coeff.*gx;  qy = -coeff.*gy;
	contrib = A .* (qx.*bx + qy.*by);                                % ne x 3: flux into each vertex of each element
	qn = accumarray(el(:), contrib(:), [nv 1]);
	q = qn(v(:));
end

function [outflow, inflow, times] = patchWaterFlux(md, T, v)
	% Water flux (m^3/s) leaving and entering the footprint of vertices v at every saved step of T.
	% The sheet flux vector is rebuilt per element from the saved fields, the way the solver forms it:
	%   SHAKTI  q = HydrologyBasalFlux * (-grad h / |grad h|)         (BasalFlux = K |grad h| per element)
	% and integrated as q.n over the perimeter of the footprint elements (elements with all vertices in v).
	% Elements that are all floating carry no flux.
	el = md.mesh.elements;  x = md.mesh.x(:);  y = md.mesh.y(:);  nv = numel(x);  ne = size(el,1);
	inV = false(nv,1);  inV(v) = true;
	elIn = all(inV(el), 2);
	if ~any(elIn), elIn = any(inV(el), 2); end                        % tiny footprint: elements touching it

	% unique edges and their one or two adjacent elements
	pairs = sort([el(:,[1 2]); el(:,[2 3]); el(:,[3 1])], 2);
	[edgeV, ~, idx] = unique(pairs, 'rows');
	eid = repmat((1:ne)', 3, 1);
	e1 = accumarray(idx, eid, [], @min);  e2 = accumarray(idx, eid, [], @max);
	e2(accumarray(idx, 1) == 1) = 0;                                  % domain boundary
	side1 = elIn(e1);  side2 = e2 > 0 & elIn(max(e2,1));
	perim = find(xor(side1, side2));                                  % exactly one adjacent element inside
	inner = e1(perim);  inner(~side1(perim)) = e2(perim(~side1(perim)));
	outer = e2(perim);  outer(~side1(perim)) = e1(perim(~side1(perim)));
	outer(outer == 0) = inner(outer == 0);                            % boundary edge: use the inside element

	% outward unit normal and length of each perimeter edge
	ex = reshape(x(edgeV(perim,:)), [], 2);  ey = reshape(y(edgeV(perim,:)), [], 2);
	L  = hypot(ex(:,2)-ex(:,1), ey(:,2)-ey(:,1));
	nx = (ey(:,2)-ey(:,1))./L;  ny = -(ex(:,2)-ex(:,1))./L;
	cx = mean(reshape(x(el(inner,:)), [], 3), 2);  cy = mean(reshape(y(el(inner,:)), [], 3), 2);
	flip = nx.*(mean(ex,2)-cx) + ny.*(mean(ey,2)-cy) < 0;
	nx(flip) = -nx(flip);  ny(flip) = -ny(flip);

	% P1 gradient per element: df/dx = sum(bx.*f(el),2), df/dy = sum(by.*f(el),2)
	x1 = x(el(:,1)); x2 = x(el(:,2)); x3 = x(el(:,3));  y1 = y(el(:,1)); y2 = y(el(:,2)); y3 = y(el(:,3));
	twoA = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);
	bx = [y2-y3, y3-y1, y1-y2] ./ twoA;  by = [x3-x2, x1-x3, x2-x1] ./ twoA;

	nt = numel(T);  outflow = nan(1,nt);  inflow = nan(1,nt);  times = [T.time];
	for i = 1:nt
		f = double(T(i).HydrologyHead(:));
		gx = sum(bx.*f(el), 2);  gy = sum(by.*f(el), 2);  gn = max(hypot(gx,gy), eps);
		coeff = double(T(i).HydrologyBasalFlux(:)) ./ gn;
		if isfield(T, 'MaskOceanLevelset') && numel(T(i).MaskOceanLevelset) == nv
			coeff(all(T(i).MaskOceanLevelset(el) <= 0, 2)) = 0;        % all-floating elements carry no flux
		end
		qx = -coeff.*gx;  qy = -coeff.*gy;
		qn = qx(outer).*nx + qy(outer).*ny;
		outflow(i) = sum(max(0, qn).*L);
		inflow(i)  = sum(max(0, -qn).*L);
	end
end

function T = fillIceFields(md, T)
	% Hydrology-only transients save no ice fields; the flux calculations read Thickness and MaskOceanLevelset
	% per output step, so fill them with the (frozen) model state where they are missing.
	if isempty(T), return, end
	if ~isfield(T, 'Thickness') || isempty(T(end).Thickness), [T.Thickness] = deal(md.geometry.thickness(:)); end
	if ~isfield(T, 'MaskOceanLevelset') || isempty(T(end).MaskOceanLevelset), [T.MaskOceanLevelset] = deal(md.mask.ocean_levelset(:)); end
end
