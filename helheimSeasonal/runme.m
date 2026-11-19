% Standalone SHAKTI simulations of Helheim Glacier forced by seasonal meltwater
% inputs to the bed, after the "Seasonal Hydrology Forcing" experiments of
% Sommers et al. (2024, https://doi.org/10.1029/2024GL109168, section 2.4.1).
%
% Unlike in the paper, the stress balance is not solved here: the ice velocity
% stays at its winter value while the subglacial hydrology evolves.

% Meltwater scenarios to run one after the other, each one starts from the
% winter base state and is saved in Helheim_<scenario>_<resolution>.mat:
%    'seasonal'       seasonal runoff below 900 m
%    'seasonal_firn'  seasonal runoff below 900 m + steady firn aquifer drainage
%    'enhanced_melt'  doubled runoff below 1500 m + doubled firn aquifer drainage
scenarios={'seasonal','seasonal_firn','enhanced_melt'};

% Number of processors and mesh resolution: change them according to your
% machine, or set np and resolution before running this script (as
% run/runTestCases.m does)
%    'low'   300 m to 3 km (7,366 elements)
%    'high'  mesh of the paper (27,913 elements)
if ~exist('np','var'), np=2; end
if ~exist('resolution','var'), resolution='low'; end

for i=1:numel(scenarios)
	scenario=scenarios{i};
	disp(['Seasonal meltwater forcing and solve (' scenario ')']);
	switch resolution
		case 'low',  md=loadmodel('HelheimWinterSpinup_300m');
		case 'high', md=loadmodel('HelheimWinterSpinup_70m');
		otherwise, error(['resolution ''' resolution ''' not supported']);
	end

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
	runoff=load('merra2_runoff.mat','ro2018');
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
	md.cluster=generic('np',np);

	md.miscellaneous.name=['Helheim_' scenario '_' resolution];
	md.verbose.solution=1;
	md=solve(md,'Transient');

	save(['Helheim_' scenario '_' resolution],'md','-v7.3');
end

