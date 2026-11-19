function runTestCases(cases,np,resolution)
%RUNTESTCASES - run SHAKTI test cases, then their diagnostics
%
%   For each test case, this runs runme.m and then diagnostics.m in the folder
%   of the test case. Results (.mat files) are saved in that folder, and
%   videos and figures in its 'diagnostics' subfolder.
%
%   Usage:
%      runTestCases(cases,np)
%      runTestCases(cases,np,resolution)
%
%      cases       'all', the name of one test case, or a cell array of names:
%                     'singleMoulin'           (about 720 hourly time steps)
%                     'helheimSeasonal'        (3 scenarios of 1 year)
%                     'thwaitesInterventions'  (spin-up, capacity and 3 experiments)
%      np          number of processors, used for all the runs
%      resolution  'low' (default) or 'high': mesh of helheimSeasonal and
%                  thwaitesInterventions (singleMoulin has one mesh only)
%
%   Examples:
%      runTestCases('all',4)
%      runTestCases('helheimSeasonal',8,'high')
%      runTestCases({'singleMoulin','thwaitesInterventions'},2,'low')
%
%   From a terminal, in this folder:
%      matlab -batch "runTestCases('all',4)"

	available={'singleMoulin','helheimSeasonal','thwaitesInterventions'};

	% Check inputs
	if nargin<2
		help runTestCases
		error('runTestCases: the test cases and the number of processors are required');
	end
	if nargin<3, resolution='low'; end
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
	root=fileparts(fileparts(mfilename('fullpath')));
	status=repmat({'done'},1,numel(cases));
	elapsed=zeros(1,numel(cases));
	for i=1:numel(cases)
		fprintf('\n=== %s (np=%d, resolution=%s) ===\n',cases{i},np,resolution);
		start=tic;
		try
			runcase(fullfile(root,cases{i}),np,resolution);
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
