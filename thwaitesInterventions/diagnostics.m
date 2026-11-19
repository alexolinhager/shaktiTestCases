% Video, stills and time series of each Thwaites run (runme.m), written in
% ./diagnostics. See ../shaktiDiagnostics.m
addpath('..');

runs={'ThwaitesSpinup','ThwaitesCapacity','Thwaites_control','Thwaites_pump','Thwaites_inject'};

% Mesh resolution of the runs ('low' or 'high'), as in runme.m
if ~exist('resolution','var'), resolution='low'; end

% Intervention sites, as in runme.m
pumpXY=[-1429020 -441873
        -1438900 -424922
        -1419150 -463010];

for i=1:numel(runs)
	name=[runs{i} '_' resolution];
	file=[name '.mat'];
	if ~exist(file,'file'), disp(['   ' file ' not found, skipped']); continue; end
	md=loadmodel(file);

	% Effective pressure is also averaged within 200 km of the grounding line
	% (where md.mask.ocean_levelset is zero)
	gl=isoline(md,md.mask.ocean_levelset,'value',0,'output','matrix');
	gl=gl(~isnan(gl(:,1)),1:2);
	[pos,distance]=dsearchn(gl,[md.mesh.x md.mesh.y]);
	region=find(distance<=200e3);
	options={'name',name,'sites',pumpXY,'region',region,'regionname','Within 200 km of the grounding line'};

	% Capacity run: the head is pinned over the sites, so what the pumps take
	% out is the net inflow into the sites computed in step 2 of runme.m
	capfile=['ThwaitesCapacity_cap_' resolution '.mat'];
	if strcmp(runs{i},'ThwaitesCapacity') && exist(capfile,'file')
		cap=load(capfile);
		options=[options {'series',struct('name','Pump output','time',cap.cap.time,'value',sum(cap.cap.Qin,1))}];
	end

	shaktiDiagnostics(md,options{:});
end
