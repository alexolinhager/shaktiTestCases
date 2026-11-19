function out=shaktiDiagnostics(md,varargin)
%SHAKTIDIAGNOSTICS - video, stills and time series of a SHAKTI run
%
%   Files written (in ./diagnostics by default):
%      <name>.mp4             effective pressure, water flux, gap height and
%                             Reynolds number at every saved time step
%      <name>_first.png       the same maps at the first saved time step
%      <name>_last.png        the same maps at the last saved time step
%      <name>_timeseries.png  water put into (or taken out of) the bed, and
%                             mean effective pressure
%
%   Only grounded ice is shown and averaged. The Reynolds number is the saved
%   one if 'HydrologyReynolds' was among the requested outputs, otherwise it
%   is the water flux divided by the kinematic viscosity of water (which is
%   how SHAKTI computes it).
%
%   Usage:
%      out=shaktiDiagnostics(md,options)
%
%   Options:
%      'name'        prefix of the files (default: md.miscellaneous.name)
%      'folder'      where the files go (default: 'diagnostics')
%      'region'      vertices over which the effective pressure is also
%                    averaged (indices or logical mask)
%      'regionname'  name of that region in the legend
%      'sites'       [x y] of points to mark on the maps (default: vertices
%                    with a moulin input, if there are 100 or fewer)
%      'series'      extra curves for the water plot: structure array with
%                    fields name, time (yr) and value (m^3/s)
%
%   Example:
%      md=loadmodel('MoulinTransient');
%      shaktiDiagnostics(md,'name','singleMoulin');

	options=pairoptions(varargin{:});
	name=getfieldvalue(options,'name',md.miscellaneous.name);
	folder=getfieldvalue(options,'folder','diagnostics');
	region=getfieldvalue(options,'region',[]);
	regionname=getfieldvalue(options,'regionname','Region');
	extra=getfieldvalue(options,'series',[]);
	if ~exist(folder,'dir'), mkdir(folder); end

	NU=1.787e-6; % kinematic viscosity of water (m^2/s), as in SHAKTI
	nv=md.mesh.numberofvertices;
	ne=md.mesh.numberofelements;
	index=md.mesh.elements;
	x=md.mesh.x(:);
	y=md.mesh.y(:);
	yts=md.constants.yts;

	% Saved time steps (a run that stopped early can end with an empty one)
	T=md.results.TransientSolution;
	keep=arrayfun(@(s) ~isempty(s.time) && numel(s.EffectivePressure)==nv && ...
		numel(s.HydrologyBasalFlux)==ne && numel(s.HydrologyGapHeight)==ne,T);
	T=T(keep);
	if isempty(T), error('shaktiDiagnostics: no hydrology results in md.results.TransientSolution'); end
	time=[T.time];
	nt=numel(time);
	N=double([T.EffectivePressure])/1e6; % MPa
	q=double([T.HydrologyBasalFlux]);    % m^2/s
	gap=double([T.HydrologyGapHeight]);  % m
	if isfield(T,'HydrologyReynolds') && all(arrayfun(@(s) numel(s.HydrologyReynolds)==ne,T))
		reynolds=double([T.HydrologyReynolds]);
	else
		reynolds=q/NU;
	end

	% Grounded ice: the elements shown, and the area that each vertex stands for
	ice=md.mask.ice_levelset(:)<0 & md.mask.ocean_levelset(:)>0;
	active=all(ice(index),2);
	if ~any(active), error('shaktiDiagnostics: no grounded ice'); end
	areas=GetAreas(index,x,y);
	weights=accumarray(reshape(index(active,:),[],1),repmat(areas(active)/3,3,1),[nv 1]);

	% Time in days for a run of up to 2 years, in years otherwise
	if time(end)-time(1)<=2
		timescale=yts/86400; timename='day';  timelabel='Time (days)';
	else
		timescale=1;         timename='year'; timelabel='Time (yr)';
	end
	t=time*timescale;

	% 1. Water put into the bed and taken out of it (m^3/s), at 500 times from
	% the start of the run (the forcing can vary faster than results are saved)
	tin=linspace(min(md.timestepping.start_time,time(1)),time(end),500);
	moulin=attimes(md.hydrology.moulin_input,nv,tin);
	englacial=attimes(md.hydrology.englacial_input,nv,tin); % m/yr
	series=struct('name',{'Moulin input','Englacial input','Pump output'}, ...
		'time',tin, ...
		'value',{sum(max(moulin,0),1), weights'*max(englacial,0)/yts, -sum(min(moulin,0),1)});
	for i=1:numel(extra)
		series(end+1)=struct('name',extra(i).name,'time',extra(i).time(:)','value',extra(i).value(:)'); %#ok<AGROW>
	end
	series=series(arrayfun(@(s) any(s.value~=0),series));

	% 2. Mean effective pressure (MPa), weighted by area
	Nmean=(weights'*N)/sum(weights);
	Nregion=[];
	if ~isempty(region)
		inregion=false(nv,1);
		inregion(region)=true;
		wr=weights.*inregion;
		if sum(wr)==0, error('shaktiDiagnostics: ''region'' holds no grounded ice'); end
		Nregion=(wr'*N)/sum(wr);
	end

	% Points to mark on the maps
	sites=getfieldvalue(options,'sites',[]);
	if ~exist(options,'sites')
		pos=find(any(moulin~=0,2));
		if numel(pos)<=100, sites=[x(pos) y(pos)]; end
	end

	% 3. Maps: same color limits for all time steps (0.1 to 99.9th percentile)
	fields=struct('name',{'Effective pressure (MPa)','Water flux (m^2 s^{-1})','Gap height (m)','Reynolds number'}, ...
		'data',{N(ice,:),q(active,:),gap(active,:),reynolds(active,:)}, ...
		'log',{false,true,true,true});
	if max(x(ice))-min(x(ice))>5e3, scale=1e-3; unit='km'; else, scale=1; unit='m'; end
	fig=figure('Visible','off','Color','w','Units','pixels','Position',[50 50 1280 960]);
	layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
	h=gobjects(4,1);
	for i=1:4
		ax=nexttile(layout);
		patch(ax,'Faces',index(~active,:),'Vertices',[x y]*scale,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
		if i==1 % on vertices
			h(i)=patch(ax,'Faces',index(active,:),'Vertices',[x y]*scale,'FaceVertexCData',N(:,1),'FaceColor','interp','EdgeColor','none');
		else    % on elements
			h(i)=patch(ax,'Faces',index(active,:),'Vertices',[x y]*scale,'FaceVertexCData',fields(i).data(:,1),'FaceColor','flat','EdgeColor','none');
		end
		if ~isempty(sites)
			line(ax,sites(:,1)*scale,sites(:,2)*scale,'LineStyle','none','Marker','o','MarkerSize',7,'Color','k','LineWidth',1);
		end
		axis(ax,'equal');
		margin=0.02*max(max(x(ice))-min(x(ice)),max(y(ice))-min(y(ice)))*scale;
		axis(ax,[min(x(ice))*scale-margin max(x(ice))*scale+margin min(y(ice))*scale-margin max(y(ice))*scale+margin]);
		set(ax,'Box','off','Layer','top','TickDir','out','FontSize',11);
		xlabel(ax,['x (' unit ')']);
		ylabel(ax,['y (' unit ')']);
		title(ax,fields(i).name,'FontWeight','normal','FontSize',13);
		values=fields(i).data(:);
		if fields(i).log
			values=values(values>0);
			if isempty(values), values=1; end
			set(ax,'ColorScale','log');
		end
		limits=prctile(values,[0.1 99.9]);
		if fields(i).log
			if limits(2)<100*limits(1) % whole decades, so that the colorbar has labels
				limits=10.^[floor(log10(limits(1))) ceil(log10(limits(2)))];
				if limits(2)<=limits(1), limits(2)=10*limits(1); end
			end
		elseif limits(2)<=limits(1)
			limits=limits(1)+[-0.5 0.5]*max(abs(limits(1)),1e-6);
		end
		clim(ax,limits);
		if ~fields(i).log && limits(1)<0 && limits(2)>0
			colormap(ax,diverging(256,limits)); % water pressure above overburden stands out
		else
			colormap(ax,sequential(256));
		end
		colorbar(ax);
	end

	% Video, with stills of the first and last time steps
	try
		video=VideoWriter(fullfile(folder,[name '.mp4']),'MPEG-4');
	catch % MPEG-4 is not available on Linux
		video=VideoWriter(fullfile(folder,[name '.avi']),'Motion JPEG AVI');
	end
	video.FrameRate=min(max(round(nt/20),2),30);
	video.Quality=90;
	open(video);
	for k=1:nt
		set(h(1),'FaceVertexCData',N(:,k));
		set(h(2),'FaceVertexCData',q(active,k));
		set(h(3),'FaceVertexCData',gap(active,k));
		set(h(4),'FaceVertexCData',reynolds(active,k));
		title(layout,sprintf('%s: %s %.2f',name,timename,t(k)),'Interpreter','none','FontSize',15);
		frame=print(fig,'-RGBImage','-r96');
		frame=frame(1:2*floor(size(frame,1)/2),1:2*floor(size(frame,2)/2),:);
		if k==1, framesize=size(frame); end
		if ~isequal(size(frame),framesize), frame=imresize(frame,framesize(1:2)); end
		writeVideo(video,frame);
		if k==1,  exportgraphics(fig,fullfile(folder,[name '_first.png']),'Resolution',150); end
		if k==nt, exportgraphics(fig,fullfile(folder,[name '_last.png']),'Resolution',150); end
	end
	close(video);
	close(fig);

	% 4. Time series: water on top (if there is any), mean effective pressure below
	colors=[42 120 214; 235 104 52; 27 175 122; 237 161 0]/255;
	fig=figure('Visible','off','Color','w','Units','pixels','Position',[50 50 900 300+300*~isempty(series)]);
	layout=tiledlayout(fig,1+~isempty(series),1,'TileSpacing','compact','Padding','compact');
	if ~isempty(series)
		ax=nexttile(layout);
		hold(ax,'on');
		for i=1:numel(series)
			plot(ax,series(i).time*timescale,series(i).value,'LineWidth',2,'Color',colors(mod(i-1,4)+1,:));
		end
		ylabel(ax,'Water (m^3 s^{-1})');
		if numel(series)>1
			legend(ax,{series.name},'Location','best','Box','off');
		else
			title(ax,series.name,'FontWeight','normal','FontSize',13);
		end
		style(ax,[tin(1) time(end)]*timescale);
	end
	ax=nexttile(layout);
	hold(ax,'on');
	plot(ax,t,Nmean,'LineWidth',2,'Color',colors(1,:));
	if ~isempty(Nregion)
		plot(ax,t,Nregion,'LineWidth',2,'Color',colors(2,:));
		legend(ax,{'All grounded ice',regionname},'Location','best','Box','off');
	end
	ylabel(ax,'Mean effective pressure (MPa)');
	xlabel(ax,timelabel);
	style(ax,[tin(1) time(end)]*timescale);
	title(layout,name,'Interpreter','none','FontSize',15);
	exportgraphics(fig,fullfile(folder,[name '_timeseries.png']),'Resolution',150);
	close(fig);

	out=struct('time',time,'series',series,'effective_pressure',Nmean,'effective_pressure_region',Nregion);
	disp(['   diagnostics of ' name ' (' num2str(nt) ' time steps) written in ' folder]);
end

function F=attimes(field,n,time) % {{{
	% Field on vertices at the given times: one column per time. A transient
	% field holds its time stamps in its last row; ISSM interpolates linearly
	% between them, and holds the first and last values outside
	if size(field,1)==n+1 && size(field,2)>1
		stamps=field(end,:);
		F=interp1(stamps(:),field(1:n,:)',min(max(time(:),stamps(1)),stamps(end)));
		F=reshape(F,numel(time),n)';
	else
		F=repmat(field(1:n,1),1,numel(time));
	end
end % }}}
function map=sequential(n) % {{{
	% One hue, from light to dark
	anchors=[234 242 251; 143 190 240; 42 120 214; 20 70 140; 8 30 70]/255;
	map=interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n));
end % }}}
function map=diverging(n,limits) % {{{
	% Two hues around a neutral gray at zero, with the same color for the same
	% distance from zero on both sides
	negative=[226 226 222; 240 160 120; 235 104 52; 150 50 15; 80 20 0]/255;
	positive=[226 226 222; 143 190 240; 42 120 214; 20 70 140; 8 30 70]/255;
	values=linspace(limits(1),limits(2),n)'/max(abs(limits));
	map=zeros(n,3);
	pos=values<0;
	map(pos,:)=interp1(linspace(0,1,size(negative,1)),negative,-values(pos));
	map(~pos,:)=interp1(linspace(0,1,size(positive,1)),positive,values(~pos));
end % }}}
function style(ax,limits) % {{{
	set(ax,'Box','off','TickDir','out','FontSize',11,'XGrid','on','YGrid','on','GridColor',[0.5 0.5 0.5],'GridAlpha',0.25);
	if limits(2)>limits(1), xlim(ax,limits); end
end % }}}
