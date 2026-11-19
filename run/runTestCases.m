function runTestCases(cases,np,rundir,resolution)
%RUNTESTCASES - run SHAKTI test cases, then their diagnostics
%
%   Nothing is written in the repository: each test case runs in its own
%   subfolder of the run directory, <rundir>/<test case>. The scripts and
%   input files of the test case are copied there, then runme.m and
%   diagnostics.m are run there. Results (.mat files) are saved in that
%   subfolder, and videos and figures in its 'diagnostics' subfolder.
%
%   Usage:
%      runTestCases(cases,np,rundir)
%      runTestCases(cases,np,rundir,resolution)
%
%      cases       'all', the name of one test case, or a cell array of names:
%                     'singleMoulin'           (about 720 hourly time steps)
%                     'helheimSeasonal'        (3 scenarios of 1 year)
%                     'thwaitesInterventions'  (spin-up, capacity and 3 experiments)
%      np          number of processors, used for all the runs
%      rundir      run directory, created if needed (outside of the repository)
%      resolution  'low' (default) or 'high': mesh of helheimSeasonal and
%                  thwaitesInterventions (singleMoulin has one mesh only)
%
%   Examples:
%      runTestCases('all',4,'~/shaktiRuns')
%      runTestCases('helheimSeasonal',8,'/scratch/shaktiRuns','high')
%      runTestCases({'singleMoulin','thwaitesInterventions'},2,'~/shaktiRuns','low')
%
%   From a terminal, in this folder:
%      matlab -batch "runTestCases('all',4,'~/shaktiRuns')"

	available={'singleMoulin','helheimSeasonal','thwaitesInterventions'};

	% Input files of each test case, copied to the run directory with runme.m
	% and diagnostics.m
	inputs.singleMoulin={'moulin.par','outline.exp'};
	inputs.helheimSeasonal={'merra2_runoff.mat','HelheimWinterSpinup_300m.mat','HelheimWinterSpinup_70m.mat'};
	inputs.thwaitesInterventions={'TIS_island.exp','Thwaites_SHAKTI_init.mat','Thwaites_SHAKTI_init_500mto10km.mat'};

	% Check inputs
	if nargin<3
		help runTestCases
		error('runTestCases: the test cases, the number of processors and the run directory are required');
	end
	if nargin<4, resolution='low'; end
	if ischar(cases) || isstring(cases)
		if strcmpi(cases,'all'), cases=available; else, cases=cellstr(cases); end
	end
	for i=1:numel(cases)
		if ~ismember(cases{i},available)
			error(['runTestCases: test case ''' cases{i} ''' not found, choose among: ' strjoin(available,', ') ', or ''all''']);
		end
	end
	if ~isnumeric(np) || ~isscalar(np) || np<1 || np~=round(np)
		error('runTestCases: np (number of processors) should be a positive integer');
	end
	if ~ismember(resolution,{'low','high'})
		error('runTestCases: resolution should be ''low'' or ''high''');
	end
	if ~(ischar(rundir) || isstring(rundir)) || strlength(rundir)==0
		error('runTestCases: rundir (run directory) should be the name of a directory');
	end

	% Run directory: absolute path, outside of the repository
	root=fileparts(fileparts(mfilename('fullpath')));
	rundir=char(rundir);
	if strncmp(rundir,'~',1), rundir=[char(java.lang.System.getProperty('user.home')) rundir(2:end)]; end
	if ~exist(rundir,'dir'), mkdir(rundir); end
	rundir=cd(cd(rundir));
	if strcmp(rundir,root) || strncmp(rundir,[root filesep],numel(root)+1)
		error(['runTestCases: the run directory should be outside of the repository (' root ')']);
	end
	addpath(root); % shaktiDiagnostics.m

	% ISSM: load its path if needed. ISSM_DIR is also needed to launch a run,
	% and is not set when MATLAB was not started from a terminal
	if ~exist('loadmodel','file')
		if isempty(getenv('ISSM_DIR'))
			error('runTestCases: ISSM is not on the MATLAB path and the ISSM_DIR environment variable is not set');
		end
		addpath([getenv('ISSM_DIR') '/src/m/dev']);
		devpath;
	end
	if isempty(getenv('ISSM_DIR')), setenv('ISSM_DIR',issmdir()); end

	% Run the test cases one after the other; if one fails, go on with the next
	status=repmat({'done'},1,numel(cases));
	elapsed=zeros(1,numel(cases));
	for i=1:numel(cases)
		fprintf('\n=== %s (np=%d, resolution=%s) in %s ===\n',cases{i},np,resolution,fullfile(rundir,cases{i}));
		start=tic;
		try
			folder=fullfile(rundir,cases{i});
			if ~exist(folder,'dir'), mkdir(folder); end
			files=[{'runme.m','diagnostics.m'} inputs.(cases{i})];
			for j=1:numel(files)
				copyfile(fullfile(root,cases{i},files{j}),folder);
			end
			runcase(folder,np,resolution);
		catch me
			status{i}='FAILED';
			disp(getReport(me));
		end
		elapsed(i)=toc(start);
	end

	% Summary
	fprintf('\n%-24s %-8s %s\n','test case','status','time (min)');
	for i=1:numel(cases)
		fprintf('%-24s %-8s %.1f\n',cases{i},status{i},elapsed(i)/60);
	end
	if any(strcmp(status,'FAILED'))
		error('runTestCases: %d of %d test cases failed (see above)',sum(strcmp(status,'FAILED')),numel(cases));
	end
end

function runcase(folder,np,resolution) % {{{
	% runme.m and diagnostics.m are scripts: they read np and resolution from
	% this workspace, and leave their variables in it
	origin=pwd;
	restore=onCleanup(@() cd(origin));
	cd(folder);
	run('runme.m');
	run('diagnostics.m');
	close all
end % }}}
