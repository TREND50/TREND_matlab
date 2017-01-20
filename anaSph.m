function recons_paras=anaSph(nrun,cid)    
% Displays amplitude as a function of distanceto emission point
% + spherical fit
% OMH 06/06/119


sharedGlobals;
scrsz = get( 0, 'ScreenSize' );

%% Plot options
%%===
labelOpts = { 'FontSize', 18 };
axisOpts  = { 'fontSize', 14, 'YScale', 'linear' };

%% Load dst
dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
dst = load(dstname);
Struct = dst.Struct;
CoincStruct = Struct.Coinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
x = [DetStruct.X];
y = [DetStruct.Y];
z = [DetStruct.Z];
DetPos = [x' y' z'];

% Coinc paras
tag = [CoincStruct.Det.Tag];
amp = [CoincStruct.Det.AmpMax];
maxraw = [CoincStruct.Det.MaxRaw];
minraw = [CoincStruct.Det.MinRaw];
sat = [maxraw==255 | minraw==0];
gain = [CoincStruct.Det.Gain];

% Recons paras
CoincId = CoincStruct.IdCoinc;
L = CoincStruct.Mult;
x0 = CoincStruct.SphRecons.X0;
y0 = CoincStruct.SphRecons.Y0;
z0 = CoincStruct.SphRecons.Z0;
rhos = CoincStruct.SphRecons.Rho;
thetas = CoincStruct.SphRecons.Theta;
phis = CoincStruct.SphRecons.Phi;
dist = Struct.Coinc.SphRecons.minDistSource;

if ~exist('cid')   
    cid = CoincId;
end

%% Loop on coincs
for i=1:length(cid)  % Scan all coincs in the list
    id = cid(i); % Coinc number
    ind = find(CoincId==id);
    if isempty(ind)
       disp(sprintf('Coinc %d not found in list of reconstructed coincs.',id))
       continue
    end
    in = find(tag(ind,:)==1);
    bsat = find(sat(ind,in)==1);
    bok = find(sat(ind,in)==0);
    calamp = [amp(ind,in).*gain(ind,in)]';
    dets = Detectors(in)';
    error = ErrorAmp*calamp;
    %
    posAnt = [x(in)' y(in)' z(in)'];
    posSource = ones(L(ind),1)*[x0(ind) y0(ind) z0(ind)];

    if isempty(ind)
          disp(sprintf('Coinc %d not found in list of reconstructed coincs.',id))
          return
    end
    disp(sprintf('Reconstructed souce position: x = %3.1f m, y = %3.1f m, z = %5.1f m',x0(ind),y0(ind),z0(ind)));
    dist = sqrt(sum((posAnt-posSource).^2,2));
    %disp 'AntennaId  RawAmp CalAmp  Dist2Source Saturation'
    %[dets amp(ind,in)' calamp dist sat(ind,in)']
    
    % Plot
    figure(9)
    title = sprintf('Amplitude vs distance to point source - Run %d Coinc %d ',nrun,id);
    set(9,'Name',title,'NumberTitle','off')
    errorbar(dist(bok),calamp(bok),error(bok),'ks','MarkerFaceColor','k')
    hold on
    plot(dist(bsat),calamp(bsat),'^r','MarkerFaceColor','r')  % Saturated events
    for j = 1:L(ind)
       text(dist(j) + 10,calamp(j),num2str(dets(j)),'FontSize',14);
    end
    axis([min(dist)*0.99 max(dist)*1.01 0 max(calamp*1.1)])
    xlabel('Distance to source [m]',labelOpts{:})
    ylabel('Signal Amplitude [LSB]',labelOpts{:})
    grid on

    % Fit
    A0 = FitInvert(dist,calamp);
    if min(dist)<1e5
        d = min(dist)*0.99:max(dist)*1.01;
        asph = A0./d;
        plot(d,asph,'r','LineWidth',2)
    else
        disp('Source too far away to plot radius.')
    end
    clear d
    clear asph
end

