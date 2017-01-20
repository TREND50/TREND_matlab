function recons_paras=AnaLDF(nrun,cid,dst)    
% Displays lateral distribution function
% OMH 07/06/2011

SharedGlobals;

%% Load dst
if ~exist('dst')
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    dst = load(dstname);
end
Struct = dst.Struct;
Struct = Dist2Source(Struct);
CoincStruct = Struct.Coinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
Type = [DetStruct.isScint];
x = [DetStruct.X];
y = [DetStruct.Y];
z = [DetStruct.Z];
DetPos = [x' y' z'];

% Coinc paras
amp = [CoincStruct.Det.AmpMax];
calampbline = Struct.Coinc.Det.CalibratedAmp1;
calamppsd = Struct.Coinc.Det.CalibratedAmp2;
evt = [CoincStruct.Det.Evt];
stat = [CoincStruct.Det.Status];
% Recons paras
CoincId = CoincStruct.IdCoinc;
L = CoincStruct.Mult;
Lant = CoincStruct.MultAnt;
Lsci = L-Lant;
thetap = CoincStruct.PlanRecons.Radio.Theta;
phip = CoincStruct.PlanRecons.Radio.Phi;

% ShowerParas
IsShower = CoincStruct.IsShower;
xCore = CoincStruct.ShowerRecons.XCore;
yCore = CoincStruct.ShowerRecons.YCore;
zCore = CoincStruct.ShowerRecons.ZCore;
AxisAmpdB = CoincStruct.ShowerRecons.AxisAmp;
AxisAmp = db2mag(AxisAmpdB);
lambda = CoincStruct.ShowerRecons.Lambda;


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
    if IsShower(ind)==0
       disp(sprintf('No shower recons for coinc %d.',id))
       continue
    end
    status = stat(ind,:);
    thisStatus = zeros(length(Detectors),10);
    thisStatus(:,1) = Detectors';
    for j = 1:length(Detectors)
      st = dec2bin(status(j));
      stv = int16(sscanf(st,'%1d'))'; %LSB last
      stv=fliplr(stv); %LSB first
      thisStatus(j,2:length(stv)+1) = stv;
    end
    in = find(thisStatus(:,2)==1);
    ant = find(Type == 0);
    sci = find(Type == 1);
    ina = intersect(in,ant);  % index for antennas in coinc
    ants = Detectors(ina)';
    scis = Detectors(sci)';
    sat = find(thisStatus(:,2)==1);
    badsig = find(thisStatus(:,3)==1);
    highrate = find(thisStatus(:,4)==1);
    
    bsat = (thisStatus(in,5)==1);
    bok = (thisStatus(in,2)==1 & thisStatus(in,5)==0);
    
    thisCalampBline = calampbline(ind,ina)'; 
    thisCalampPsd = calamppsd(ind,ina)';
    if max(thisCalampPsd)>0        
        isPSD = 1;
    else
        isPSD = 0;  % No psd measurement
    end
    error = ErrorAmp*thisCalampBline;
    errorpsd = ErrorAmp*thisCalampPsd;
    
    % Determine detector distance to shower axis
    posDet = [x' y' z'];
    posAnt = [x(ina)' y(ina)' z(ina)'];
    posCore = (ones(length(Detectors),1)*[xCore(ind) yCore(ind) zCore(ind)]);
    posCore_ant = ones(Lant(ind),1)*[xCore(ind) yCore(ind) zCore(ind)];
    dX =  posDet - posCore;
    %dX =  posAnt - posCore_ant;
    cp = cosd(phip(ind)); 
    sp = sind(phip(ind));
    ct = cosd(thetap(ind)); 
    st = sind(thetap(ind));
    u  = [ -sp*st, cp*st, ct ];  %% Warning!!!! This is different from usual since here x=WE & y=SN
    U  = ones( length(Detectors), 1 )*u;
    
    dist = sqrt( sum( dX.^2, 2 ) - sum( dX.*U, 2 ).^2 ); %Distance from antennas to shower axis     
    dist_sci = dist(sci);
    dist_ant = dist(ina);    
    
    if isempty(ind)
          disp(sprintf('Coinc %d not found in list of reconstructed coincs.',id))
          return
    end
    disp(sprintf('Reconstructed core position: x = %3.1f m, y = %3.1f m, z = %5.1f m',xCore(ind),yCore(ind),zCore(ind)));
    disp(sprintf('Shower amplitude on axis = %3.1f LSB, Lamba = %3.1f m',AxisAmp(ind),lambda(ind)))
    disp(sprintf('ScintId Dist to shower axis  Amplitude '))
    for j=1:length(scis)
        disp(sprintf('   %d        %3.1f m         %3.2f V',scis(j),dist_sci(j),amp(ind,sci(j))))    
    end
    
    % Plot
    figure(8)
    title = sprintf('Amplitude Lateral Distribution Function - Run %d Coinc %d ',nrun,id);
    set(8,'Name',title,'NumberTitle','off')
    errorbar(dist_ant(bok),thisCalampBline(bok),error(bok),'ks','MarkerFaceColor','w')
    plot(dist_ant(bsat),thisCalampBline(bsat),'^r','MarkerFaceColor','r')  % Saturated events
    hold on
    if isPSD == 1
        errorbar(dist_ant(bok),thisCalampPsd(bok),errorpsd(bok),'ks','MarkerFaceColor','k')  
        ylabel('Signal Amplitude [V]',labelOpts{:})
        for j = 1:Lant(ind)
            text(dist_ant(j) + 10,thisCalampPsd(j),num2str(ants(j)),'FontSize',14);
        end
    else
        ylabel('Signal Amplitude [A.U.]',labelOpts{:})
        for j = 1:Lant(ind)
            text(dist_ant(j) + 10,thisCalampBline(j),num2str(ants(j)),'FontSize',14);
        end        
    end
    %FitExp(dist(bok),thisCalampBline(bok),[0.1 1/800])
    xlabel('Distance to shower axis [m]',labelOpts{:})
    grid on

    % Fit
    if max(dist)<1e4
        d = 0:max(dist)*1.01;
        amp_th = AxisAmp(ind).*exp(-d./lambda(ind));
        plot(d,amp_th,'k','LineWidth',2)
        minx = min(dist_ant)*0.9;
        maxx = max(dist_ant)*1.1;
        miny = AxisAmp(ind).*exp(-minx./lambda(ind));
        maxy = AxisAmp(ind).*exp(-maxx./lambda(ind))*2;
        if isPSD
            text(min(dist),max(thisCalampPsd)*0.2,sprintf('lambda=%3.1f m',lambda(ind)),'FontSize',14) 
            axis([minx maxx 0 max([maxy,max(thisCalampPsd)])*1.1])
        else
            text(min(dist),max(thisCalampBline)*0.2,sprintf('lambda=%3.1f m',lambda(ind)),'FontSize',14) 
            axis([minx maxx 0 max([maxy,max(thisCalampBline)])*1.1])
        end
        clear amp_th d
    end
end

