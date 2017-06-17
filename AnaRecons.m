function list = AnaRecons(nrun,reft,dstf)
% Antennas coincidence reconstruction results
% OMH  04/01/2011


SharedGlobals;
mtyp= ['+bla';'.blu';'dmag';'pred';'hgre';'oyel'];
htyp = ['sb';'sk';'sr'];
linf = 4;
lsup = 14;
step = 2;
    
%% Load dst
if ~exist('dstf')
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    %dstname = '/media/dst/sel/dst256104_1.mat';
    dst = load(dstname);
else
    dst = dstf;
end
dst.Struct = Dist2Source(dst.Struct);
%dst.Struct = CleanFilt(dst.Struct);

ncoincs = dst.Struct.Setup.TotalCoinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
DetectorType = [dst.Struct.Setup.Det.isScint];
evt = [dst.Struct.Coinc.Det.Evt];
nDets = length(Detectors);
nScints = sum(DetectorType);
nAnts = nDets-nScints;
Events = [DetStruct.Evt];
CoincStruct = dst.Struct.Coinc;
cordelays = 0;
if cordelays == 1
    PlanStruct = CoincStruct.DelayCorrRecons.PlanRecons;
    SphStruct = CoincStruct.DelayCorrRecons.SphRecons;
else
    PlanStruct = CoincStruct.PlanRecons.Radio;
    SphStruct = CoincStruct.SphRecons;
end
PlanStructIni = CoincStruct.PlanRecons.Radio;
thetaip = PlanStructIni.Theta;
phiip = PlanStructIni.Phi;
SphStructIni = CoincStruct.SphRecons;
thetais = SphStructIni.Theta;
phiis = SphStructIni.Phi;

RunSetup = dst.Struct.Setup;
DetPosX = [RunSetup.Det.X];
DetPosY = [RunSetup.Det.Y];
DetPosZ = [RunSetup.Det.Z];

mult = CoincStruct.Mult;
time = max(CoincStruct.Det.UnixTime,[],2)/60;  % minutes
date_start = min(time);
time = time-date_start;
idp = CoincStruct.IdCoinc;
thetap = PlanStruct.Theta;
phip = PlanStruct.Phi;
thetas = SphStruct.Theta;
phis = SphStruct.Phi;
rhos = SphStruct.Rho;
xs = SphStruct.X0;
ys = SphStruct.Y0;
zs = SphStruct.Z0;
r = SphStructIni.minDistSource;
rcorr = SphStruct.minDistSource;

d = thetas-thetap;
dp = phis-phip;

%bary = CoincStruct.DistBary;
chi2 = SphStruct.Chi2Delay;
slope = SphStruct.SlopeDelay;

if exist('reft') & reft>0
    disp 'Selecting events in time window.'
    sel = find(mult>0 & time>reft-2 & time<reft+2 );
else
    sel = find(mult>=4);
    %sel = find(mult>=6 & time<18 & time>14 & thetas'<85 & thetap'<85); %Plane track R3005
    %sel = find(mult>=6 & time<18 & time>14 & thetas'<85 & thetap'<85  & abs(d')<5 & chi2'<50 & phis'<150 & phis'>0); %Plane track R3005
    %sel = find(mult>5 & time<995 & time>985 & thetas'<85 & thetap'<85  & chi2'<20 & (phis'<180 | phis'>350)); %Plane track R3004
    %sel = find(time'<7 & time'>0 & phis<150 & thetas<80 & thetap<80 & mult'>=6 & abs(d)<5 & r'>1500 & rcorr'>1500 & thetais<80 & phiis<150 & abs(thetais-thetaip)<5); %Track R256104
    %sel = find(mult'>=6 & time'<6.5 & time'>1.3 & phis<150 & thetas<78 & thetap<78 & abs(d)<5 & r'>2000  & chi2<50 & abs(slope-1)<0.1); %Track R256104
    %sel = find(mult'>=6 & time'<375 & time'>355 & phis<100 & thetas<78 & thetap<78 & abs(d)<5 & r'>2000  & chi2<50 & abs(slope-1)<0.1); %Track R3036
    %sel = find(mult'>=6 & time'<941 & time'>935 & phis<200 & thetas<78 & thetap<78 & abs(d)<5 & r'>2000  & chi2<50 & abs(slope-1)<0.1); %Track R3039  
    %sel = find(mult>=6 & chi2'<50 & abs(slope'-1)<0.1 & xs'>1400 & xs'<1600 & ys'<0 & ys'>-400);  % transfo HV
    %sel = find(mult>=6 & chi2'<50 & abs(slope'-1)<0.1 & xs'>1600 & xs'<1850 & ys'<700 & ys'>300 & time>1160);  % trace R3041
end

%PlotDelays(nrun,idp(sel))
%AnaCoinc(nrun,idp(sel))
%figure(123)
%hist(Detectors(sel))

date_start = max(dst.Struct.Setup.RunTimeStart);
[year, month, day, hour, minute, second]=unixsecs2date(date_start);
disp(sprintf('*** Run %d - %02d/%02d/%d', nrun, day, month, year))
disp(sprintf('%d antennas and %d scintillators',nAnts,nScints))
disp(sprintf('%d coincs, %d valid (%3.1f pc). Mean multiplicity=%3.1f, Max mult=%d.',ncoincs,length(sel),length(sel)/ncoincs*100,mean(mult(sel)),max(mult(sel))))
if length(sel)==0
    return
end
mult = mult(sel);
time = time(sel);
tmin = min(time);
tmax = max(time);
thetap = thetap(sel);
phip = phip(sel);
thetas = thetas(sel);
phis = phis(sel);
rhos = rhos(sel);
xs = xs(sel);
ys = ys(sel);
zs = zs(sel);
idp = idp(sel);
chi2 = chi2(sel);
slope = slope(sel);
r = r(sel);
sky = (thetas<70);
list = idp;
%list = idp(sky);
temp = zeros(length(sel),length(Detectors)); 
temp =  evt(sel,:);
evt = temp;
format long
time(1)+date_start
pause
% %% Selection to file
selection = [];
for i = 1:length(sel)
    in = find(evt(i,:)>0);
    l = length(in);
    selection(end+1:end+l,1:3) = [idp(i)*ones(l,1) Detectors(in)'  evt(i,in)'];
end
filename = sprintf('R%d_sel.mat',nrun);
save(filename,'selection')
for i = 1:length(Detectors)
  thisDet = find(selection(:,2)==Detectors(i));
  if length(thisDet)>0
    disp(sprintf('%d events saved for detector %d. 1st is %d, last is %d.',length(thisDet),Detectors(i),selection(thisDet(1),3),selection(thisDet(end),3)))
  end
end

%% Recons Analysis
figure(1)
set(1,'Name',sprintf('R%d - Multiplicity',nrun),'NumberTitle', 'off')
subplot(2,1,1)
hist(mult)
xlabel('Multiplicity', labelOpts{:})
subplot(2,1,2)
plot(chi2,slope,'+k')
xlabel('Chi2', labelOpts{:})
ylabel('Slope', labelOpts{:})
statsChi2 = datastats(chi2)
statsSlope = datastats(slope)
% Chi2 as a function of antennas in event
for i = 1:length(Detectors)
   indant = find(Detectors==Detectors(i));
   in = find(evt(:,indant)>0);
%    if Detectors(i)==132
%        thetas(in) = 0;
%        phis(in) = 0;
%    end
   mchi2(i) = mean(chi2(in));
   albl{ i } = num2str( Detectors(i) );
end
figure(8)
plot(1:length(Detectors),mchi2,'sk','MarkerFaceColor','k')
set( gca, 'XTick', 1:length(Detectors), 'XTickLabel', albl );
grid on
xlabel('Detector Id')
ylabel('Mean Chi2')
%pause

%% Azimuth vs time
figure(222)
set(222,'Name',sprintf('R%d - Azimuth vs Time',nrun),'NumberTitle', 'off','Position',[1 scrsz(2) scrsz(3)*0.95 scrsz(4)/1.2]);
subplot(2,1,1)
% Event rate
for i = 0:floor(tmax)
    eventrate(i+1) = length(find(time>i & time<i+1));
end
eventrate = eventrate/60;
plot(0:tmax,eventrate,'k','LineWidth',2)
xlabel('Time [mn]', labelOpts{:})
ylabel('Coinc rate [Hz]', labelOpts{:})
xlim([tmin tmax])
grid on
subplot(2,1,2)
ic = 1;
plot(time(sky),phis(sky),'sm','MarkerSize',6);
hold on
for lcut=linf:step:lsup
    selmult = find( mult >= lcut & mult<lcut+step);
    plot(time(selmult),phis(selmult),mtyp(ic,:),'MarkerSize',6,'LineWidth',2);
    ic=ic+1;
end
if exist('reft')
    line([reft reft],[0 360])
end
xlim([tmin tmax])
ylim([0 360])
grid on
%list = eventrate;

%% Skyplots
skyplot = 3;
prepareSkyplot(skyplot);
tit=sprintf('R%d Skyplot - Plane reconstruction',nrun);
set(skyplot,'Name',tit,'NumberTitle','off');
ic=1;
for lcut=linf:step:lsup
    selmult = find( mult >= lcut & mult<lcut+step);
    polar( phip(selmult )*DEG2RAD(1), thetap(selmult), mtyp(ic,:) );  
    hold on
    ic = ic+1;
end
%polar( phip(sky)*DEG2RAD(1), thetap(sky), 'sm' );  

skyplot = 4;
prepareSkyplot(skyplot);
tit=sprintf('R%d Skyplot - Spherical reconstruction',nrun);
set(skyplot,'Name',tit,'NumberTitle','off');
ic=1;
for lcut=linf:step:lsup
    selmult = find( mult >= lcut & mult<lcut+step);
    polar( phis(selmult )*DEG2RAD(1), thetas(selmult), mtyp(ic,:) );  
    hold on
    ic = ic+1;
end
%polar( phis(sky)*DEG2RAD(1), thetas(sky), 'sm');  

%% Details plots
close = find(r<500);
far = find(r>500);
figure(30)
set(30,'Name',sprintf('R%d - distance to source',nrun),'NumberTitle','off')
subplot(2,3,1)
hist(r(close),100)
xlabel('Distance to source [m]',labelOpts{:})
title('Spherical recons - Close sources',labelOpts{:})
subplot(2,3,2)
hist(thetas(close),100)
xlim([0 90])
xlabel('Zenith [deg]',labelOpts{:})
title('Spherical recons - Close sources',labelOpts{:})
subplot(2,3,3)
hist(phis(close),360)
xlim([0 360])
xlabel('Azimuth [deg]',labelOpts{:})
title('Spherical recons - Close sources',labelOpts{:})
%
subplot(2,3,4)
hist(log(r(far))/log(10),100)
xlabel('Distance to source (log[m]/log(10))',labelOpts{:})
title('Spherical recons - Distant sources',labelOpts{:})
subplot(2,3,5)
hist(thetas(far),100)
xlim([0 90])
xlabel('Zenith [deg]',labelOpts{:})
title('Spherical recons - Distant sources',labelOpts{:})
subplot(2,3,6)
hist(phis(far),360)
xlim([0 360])
xlabel('Azimuth [deg]',labelOpts{:})
title('Spherical recons - Distant sources',labelOpts{:})   

figure(31)
set(31,'Name',sprintf('R%d - plane reconstruction',nrun),'NumberTitle','off')
subplot(1,2,1)
hist(thetap,100)
xlim([0 90])
xlabel('Zenith [deg]',labelOpts{:})
subplot(1,2,2)
hist(phip,360)
xlim([0 360])
xlabel('Azimuth [deg]',labelOpts{:})
    
%% Position plots
figure(6)
set(6,'Name',sprintf('R%d - Reconstructed source position',nrun),'NumberTitle', 'off','Position',[1 scrsz(2) scrsz(3)*0.95 scrsz(4)/2.2]);
for i = 1:length( DetPosX )        
    text( DetPosX(i)+20, DetPosY(i), num2str( Detectors(i) ), 'FontSize', 8, 'FontWeight', 'bold','Color','k' );
end
ic=1;
for lcut=linf:step:lsup
    selmult = find( mult >= lcut & mult<lcut+step);
    plot(xs(selmult),ys(selmult),mtyp(ic,:),'MarkerSize',6,'LineWidth',2)
    ic = ic+1;
    hold on
end
plot( DetPosX(DetectorType==0), DetPosY(DetectorType==0), '^b', 'MarkerSize', 8, 'MarkerFaceColor','w','LineWidth',2 );  %Ground view  
plot( DetPosX(DetectorType==1), DetPosY(DetectorType==1), 'rs', 'MarkerSize', 8, 'MarkerFaceColor','r' );  %Ground view  
axis([-1000 5000 -2000 7000])
xlabel('W-E [m]',labelOpts{:})
ylabel('S-N [m]',labelOpts{:})
grid on
crossPointEnvir

figure(7)
set(7,'Name',sprintf('R%d - Reconstructed source position histos',nrun),'NumberTitle', 'off');
subplot(2,1,1)
hist(xs,100)
datastats(xs)
subplot(2,1,2)
hist(ys,100)
datastats(ys)
