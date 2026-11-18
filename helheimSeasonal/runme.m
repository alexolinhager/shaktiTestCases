% Standalone SHAKTI simulations of Helheim Glacier forced by seasonal meltwater
% inputs to the bed, after the "Seasonal Hydrology Forcing" experiments of
% Sommers et al. (2024, https://doi.org/10.1029/2024GL109168, section 2.4.1).
%
% Unlike in the paper, the stress balance is not solved here: the ice velocity
% stays at its winter value while the subglacial hydrology evolves.
%
% Step 1: Winter base state (initial conditions)
% Step 2: Seasonal meltwater forcing and solve, for each scenario
%
% Input files, from the archive of the paper (https://doi.org/10.5281/zenodo.11373409):
%    Helheim_SHAKTI_N_2_1yr_H10.mat   winter spin-up of coupled SHAKTI-ISSM
%    merra2_runoff.mat                14-day smoothed MERRA-2 runoff

steps=[1:2];

% Meltwater scenarios to run one after the other, each one starts from the
% winter base state and is saved in Helheim_<scenario>.mat:
%    'seasonal'       seasonal runoff below 900 m
%    'seasonal_firn'  seasonal runoff below 900 m + steady firn aquifer drainage
%    'enhanced_melt'  doubled runoff below 1500 m + doubled firn aquifer drainage
scenarios={'seasonal','seasonal_firn','enhanced_melt'};
path_to_datasets = '/Users/ahager/Documents/shaktiTestCases_datasets';

if any(steps==1)
	disp('	Step 1: Winter base state');

	% Load the winter spin-up and keep its final state
	md=loadmodel([path_to_datasets '/Helheim_SHAKTI_N_2_1yr_H10']);
	winter=md.results.TransientSolution(end);

	% Starting conditions of the subglacial system
	md.hydrology.head=winter.HydrologyHead;
	md.hydrology.gap_height=winter.HydrologyGapHeight;
	md.hydrology.reynolds=winter.HydrologyBasalFlux./1.787e-6;
	md.hydrology.reynolds(md.hydrology.reynolds==0)=1;

	% Ice sliding velocity (m/yr), held at the winter state throughout the simulation
	md.initialization.vx=winter.Vx;
	md.initialization.vy=winter.Vy;
	md.initialization.vel=winter.Vel;

	% Effective pressure at the start of the simulation, afterwards the friction
	% law uses the effective pressure computed by SHAKTI (md.friction.coupling=4)
	md.friction.effective_pressure=winter.EffectivePressure;

	% Convergence criteria of the SHAKTI nonlinear solver: the spin-up was run
	% when SHAKTI still used the stress balance criteria, so use the same ones
	md.hydrology.maxiter=md.stressbalance.maxiter;
	md.hydrology.restol=md.stressbalance.restol;
	md.hydrology.reltol=md.stressbalance.reltol;
	md.hydrology.abstol=md.stressbalance.abstol;

	% Standalone SHAKTI: only solve for hydrology
	md.transient=deactivateall(md.transient);
	md.transient.ishydrology=1;

	% Remove the output of the spin-up
	md.results=struct();

	save HelheimWinter md
end

if any(steps==2)
	for i=1:numel(scenarios)
		scenario=scenarios{i};
		disp(['	Step 2: Seasonal meltwater forcing and solve (' scenario ')']);
		md=loadmodel('HelheimWinter');

		switch scenario
			case 'seasonal',      melt_factor=1; melt_elevation=900;  firn_factor=0;
			case 'seasonal_firn', melt_factor=1; melt_elevation=900;  firn_factor=1;
			case 'enhanced_melt', melt_factor=2; melt_elevation=1500; firn_factor=2;
			otherwise, error(['scenario ''' scenario ''' not supported']);
		end

		% Define the time stepping scheme: run for 1 year with a time step of 6 hr,
		% and save the model output every 10 days
		md.timestepping.start_time=0;
		md.timestepping.time_step=6*3600/md.constants.yts; % 6 hr Time step (in years)
		md.timestepping.final_time=1;
		md.settings.output_frequency=4*10; % output every 10 days
		time=md.timestepping.time_step*(0:round(md.timestepping.final_time/md.timestepping.time_step));

		% Meltwater input rate (m/yr) from hourly 2018 MERRA-2 runoff near the
		% terminus, starting on January 1st and going back to zero at the end of the year
		runoff=load([path_to_datasets '/merra2_runoff.mat'],'ro2018');
		runoff=[runoff.ro2018(:); 0]';
		runoff=melt_factor*interp1((0:numel(runoff)-1)/(numel(runoff)-1),runoff,time);

		% Distributed input to the bed (m/yr) at low elevation, where the ice is
		% heavily crevassed
		pos=find(md.geometry.surface<=melt_elevation);
		md.hydrology.englacial_input=zeros(md.mesh.numberofvertices+1,numel(time));
		md.hydrology.englacial_input(end,:)=time;
		md.hydrology.englacial_input(pos,:)=repmat(runoff,numel(pos),1);

		% Steady drainage of the firn aquifer (50e6 m^3/yr) through crevasses at
		% 1500 m elevation, divided evenly between point inputs (m^3/s)
		pos=find(md.geometry.surface>=1500 & md.geometry.surface<=1515);
		md.hydrology.moulin_input=zeros(md.mesh.numberofvertices,1);
		md.hydrology.moulin_input(pos)=firn_factor*50e6/md.constants.yts/numel(pos);

		% Specify that you want to run the model on your current computer
		% Change the number of processors according to your machine (here np=8)
		md.cluster=generic('np',12);

		md.miscellaneous.name=['Helheim_' scenario];
		md.verbose.solution=1;
		md=solve(md,'Transient');

		save(['Helheim_' scenario],'md','-v7.3');
	end
end
