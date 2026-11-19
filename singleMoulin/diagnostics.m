% Video, stills and time series of the single moulin run (runme.m), written in
% ./diagnostics. See ../shaktiDiagnostics.m
addpath('..');

md=loadmodel('MoulinTransient');

% The moulin, at x=500, y=500
[a,pos]=min(sqrt((md.mesh.x-500).^2+(md.mesh.y-500).^2));

% Effective pressure is also averaged downstream of the moulin, between it and
% the outflow at the left edge of the domain
region=find(md.mesh.x<=md.mesh.x(pos));

shaktiDiagnostics(md,'name','singleMoulin','sites',[md.mesh.x(pos) md.mesh.y(pos)],...
	'region',region,'regionname','Downstream of the moulin');
