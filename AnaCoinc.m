function []=AnaCoinc(nrun,cid)
% Display all valuable info for coinc 'cid' in run 'irun'
% OMH 03/06/2011

SharedGlobals;

%% Load dst
%dstname = [HYB_PATH sprintf(dst_filename,nrun,1)];   % Hybrid dst
dstname = [DST_PATH sprintf(dst_filename,nrun,1)];   % Standard dst
%dstname = [DST_PATH '../hybrid/' sprintf(dsthyb_filename,nrun,1)];  % Hybrid dst
dst = load(dstname);
Struct = dst.Struct;
Struct = Dist2Source(Struct);
CoincStruct = Struct.Coinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
DetectorType = [DetStruct.isScint];
Scints = Detectors(DetectorType==1);
Ants = Detectors(DetectorType==0);
x = [DetStruct.X];
y = [DetStruct.Y];
z = [DetStruct.Z];
DetPos = [x' y' z'];

% Coinc paras
stat = [CoincStruct.Det.Status];
amp = [CoincStruct.Det.AmpMax];
calampbline = Struct.Coinc.Det.CalibratedAmp1;
calamppsd = Struct.Coinc.Det.CalibratedAmp2;
time = max(CoincStruct.Det.UnixTime,[],2)/60;  % minutes
time = time-min(time);
duration = max(time);

% starts = [Struct.Setup.InfosRun.TimeStart];
% stops = [Struct.Setup.InfosRun.TimeStop];
% stats = datastats(starts);
% good = find(abs(starts-stats.median)<3*stats.std);
% if length(good)<10
%     disp('Error in duration computing!')
%     return
% end
% start = min(starts(good));
% stop = max(stops(good));
% secs = unixs-start;
% duration = (stop-start)/60; %mins

% Recons paras
CoincId = CoincStruct.IdCoinc;
IsShower = CoincStruct.IsShower;
L = CoincStruct.Mult;
Lant = CoincStruct.MultAnt;
Lsci = L-Lant;
x0 = CoincStruct.SphRecons.X0;
y0 = CoincStruct.SphRecons.Y0;
AllDist = Struct.Coinc.SphRecons.DistSource;
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
    coinc_min = time(ind);
    status = stat(ind,:);
    thisStatus = zeros(length(Detectors),10);
    thisStatus(:,1) = Detectors';
%     tag = zeros(1,length(Detectors));
%     filtered = zeros(1,length(Detectors));
%     highrate = zeros(1,length(Detectors));
%     sat = zeros(1,length(Detectors));
%     calib = zeros(1,length(Detectors));
%     weak = zeros(1,length(Detectors));
%     noisy = zeros(1,length(Detectors));
%     envir = zeros(1,length(Detectors));
%     jump = zeros(1,length(Detectors));
    
    for j = 1:length(Detectors)
      st = dec2bin(status(j));
      stv = int16(sscanf(st,'%1d'))'; %LSB last
      stv=fliplr(stv); %LSB first
      thisStatus(j,2:length(stv)+1) = stv;
    end
    thisStatus
    in = find(thisStatus(:,2)==1);
    ramp = calamppsd(ind,in)./max(calamppsd(ind,in));
    badsig = (thisStatus(:,3)==1);
    highrate = (thisStatus(:,4)==1);
    bsat = (thisStatus(:,5)==1);
    psdin = (thisStatus(:,6)==1);
    weak = (thisStatus(:,7)==1);
    envir = (thisStatus(:,8)==1);
    noisy = (thisStatus(:,9)==1);
    jump = (thisStatus(:,10)==1);

    dets = Detectors(in);
    dpos = DetPos(in,:);
    disp(sprintf('\nCoinc %d (L=%d, Lsci =%d) at %3.1f min / %3.1f mins', id, L(ind),Lsci(ind),coinc_min,duration));  
    disp 'AntennaId   RawAmp     CalAmp      CalAmpPSD     Status'
    for j=1:L(ind)
        disp(sprintf('   %d     %3.2f V    %3.2f A.U.     %3.2f muV      %d',dets(j),amp(ind,in(j)),calampbline(in(j)),calamppsd(in(j)),status(in(j))))   
    end
    
    if isempty(ind)
          disp(sprintf('Coinc %d not found in list of reconstructed coincs.',id))
          return
    end
    
    %% Delays
    if Lant(ind)>=3  % Antennas are in
        if Lsci>1  % Hybrid
          PlotDelaysHyb(nrun,id,dst);
        else
          PlotDelays(nrun,id,dst);
        end
    else
        if Lsci(ind)==3
          disp 'This is a pure scintillator event... Skipping it.'
          continue
        else
            disp(sprintf('Error! Lsci = %d & Lant = %d',Lsci(ind), Lant(ind)))
            return
        end
    end
    
    %% Envir
    AnaRecons(nrun,coinc_min,dst);   % Skip this when using hybrid dst...
    
    %% Ground plot
    ifig = Layout(nrun,dst);
    figure(ifig)
    hold on
    % Status
    plot(DetPos(bsat,1),DetPos(bsat,2),'r^','MarkerSize',16,'MarkerFaceColor','r')
    plot(DetPos(badsig,1),DetPos(badsig,2),'kv','MarkerSize',8,'MarkerFaceColor','w')
    plot(DetPos(highrate,1),DetPos(highrate,2),'kx','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','w')
    plot(DetPos(highrate,1),DetPos(highrate,2),'kx','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','w')
    %
    plot(DetPos(psdin,1),DetPos(psdin,2),'^k','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','g')
    plot(DetPos(weak,1),DetPos(weak,2),'^k','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','k')
    plot(DetPos(envir,1),DetPos(envir,2),'^k','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r')  % Bad envir
    plot(DetPos(noisy,1),DetPos(noisy,2),'ok','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','r')  % fluctuation in bline close in time
    plot(DetPos(jump,1),DetPos(jump,2),'ok','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','b')  % jumps in this run
    %
    plot(x0(ind),y0(ind),'hg','MarkerFaceColor','g')
    [v indm] = min(AllDist(ind,find(AllDist(ind,in)>0)));
    distGround = norm([x0(ind)-dpos(indm,1) y0(ind)-dpos(indm,2)]);  % Source sphere intersection with ground at time of 1st trigger
    if dist(ind)<1e5
        myCircle(x0(ind),y0(ind),distGround,'g') 
    else
        disp 'Source too far away for radius display.'
    end
    
    for j = 1:length(dets)
        myCircle(dpos(j,1),dpos(j,2),ramp(:,j)*100)
    end
   
    %% Amplitude vs distance
    %AnaSph(nrun,id,dst)
    if IsShower(ind)  % Shower reconstruction was performed
        AnaLDF(nrun,id,dst)
        xCore = CoincStruct.ShowerRecons.XCore;
        yCore = CoincStruct.ShowerRecons.YCore;
        figure(ifig)
        plot(xCore(ind),yCore(ind),'gh','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','g')
        plotCandidate(nrun,id,dst)
    else
        disp 'This coinc is not a shower candidate.'
    end
    %     
    %%
    if length(cid)>1
        pause
        close all
    end

end    