% Video, stills and time series of each Helheim scenario (runme.m), written in
% ./diagnostics. See ../shaktiDiagnostics.m
addpath('..');

scenarios={'seasonal','seasonal_firn','enhanced_melt'};

% Mesh resolution of the runs ('low' or 'high'), as in runme.m
if ~exist('resolution','var'), resolution='low'; end

for i=1:numel(scenarios)
	name=['Helheim_' scenarios{i} '_' resolution];
	file=[name '.mat'];
	if ~exist(file,'file'), disp(['   ' file ' not found, skipped']); continue; end
	md=loadmodel(file);

	% Effective pressure is also averaged below 900 m elevation, where the
	% seasonal runoff reaches the bed in all scenarios
	region=find(md.geometry.surface<=900);

	shaktiDiagnostics(md,'name',name,...
		'region',region,'regionname','Below 900 m elevation');
end
