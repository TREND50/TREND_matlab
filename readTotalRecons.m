function [] = readTotalRecons(pstart,pstop,rstart,rstop)
% OMH 05/06/2017

SharedGlobals;

periods(1,:) = [2538 2585];  % 
periods(2,:) = [2685 2890];  %
periods(3,:) = [3000 3086];  % 
periods(4,:) = [3157 3256];  % Run -3202-4
periods(5,:) = [3336 3371];  % 
periods(6,:) = [3562 3733];  % New DAQ 
periods(7,:) = [3835 3999];  % 
periods(8,:) = [4183 4389];  % 
periods(9,:) = [4444 5014];  % NS polar  
periods(10,:) = [5070 5913]; % NS polar + DAQ upgrade  To Be Done with dst122013

if pstart == 0 & pstop == 0
  pstart = find(periods(:,1)==rstart);
  pstop = find(periods(:,2)==rstop);
end
if ~exist('pstart')
    pstart=6;
end
if ~exist('pstop')
    pstop=pstart;
end

dayToSecs = 3600*24;
dur = [22.8 42.2 44.1 49.6 23.2 79.2 26.2 28.6 58.4 63.1];   % Taken from python monitoring tool for period >=6 (/sps/hep/TREND/readSummaryProd)
cands = [44 94 87 44 54 205 24 22 33 47];
candRate = cands(pstart:pstop)./dur(pstart:pstop);
eCandRate = sqrt(cands(pstart:pstop))./dur(pstart:pstop);

col = {'k','b','g','r','m','c','y','b--','k','b'};
phitot = zeros(1,360);
thetatot = zeros(1,100);
phitotew = zeros(1,360);
thetatotew = zeros(1,100);
phitotsn = zeros(1,360);
thetatotsn = zeros(1,100);
xytot = zeros(250,250);
thetaphitot = zeros(100,360);

ew = linspace(-1000,-1000+20*250,250);
sn = linspace(-1500,-1500+20*250,250);
radtot = zeros(1,100);

nrawtot = 0;
nconscoinctot = 0;
nbadtot = 0;
rawRate = zeros(1,pstop-pstart+1);
selRate = zeros(1,pstop-pstart+1);
eRawRate = zeros(1,pstop-pstart+1);
eSelRate = zeros(1,pstop-pstart+1);


for i=pstart:pstop
    %totdst = sprintf('TotalRecons_Period%d_102014.mat',i);
    totdst = sprintf('TotalRecons_Period%d.mat',i);
    %totdst = sprintf('TotalRecons_R%dR%d_102014.mat',rstart,rstop);
    disp(sprintf('Loading %s',totdst))
    totdst = [CAND_PATH totdst];
    t = load(totdst);
    runs = t.Runs;
    durtot= sum(dur(pstart:pstop));
    [uruns iruns ] = unique(runs);   % All Setup info is DUPLICATED from one DST to the other in the same run, except for total coinc
    nraw = sum(t.CoincRaw(iruns)); % duplicated info
    nrawtot = nrawtot+nraw;
    nconscoinc = sum(t.CoincConsCoinc(iruns));  % duplicated info
    nconscoinctot = nconscoinctot+nconscoinc;
    nbad = sum(t.CoincBadPulses);  % non duplicated info
    nbadtot = nbadtot+nbad;
    disp(sprintf('Period %d: %d raw coincs, %d after ConsCoinc; %d after BadPulses',i,nraw,nconscoinc,nbad))
    %
    [nraw,nconscoinc,nbad]'
    %
    rawRate(i) = nraw./(dur(i)*dayToSecs);
    eRawRate(i) = sqrt(nraw)./(dur(i)*dayToSecs);
    selRate(i) = nbad./(dur(i)*dayToSecs);
    eSelRate(i) = sqrt(nbad)./(dur(i)*dayToSecs);
    %
    thetahist = t.ThetaSTot;
    thetatot = thetatot + thetahist;
    if i<7
        thetatotew = thetatotew+thetahist;
    else
        thetatotsn = thetatotsn+thetahist;
    end    
%     figure(2)
%     subplot(2,1,1)
%     semilogy(thetahist,col{i},'LineWidth',2)
%     %hold on
%     grid on
%     xlim([0 90])
%     ylabel('dN/d\theta [deg^{-1}]', labelOpts{:})
%     subplot(2,1,2)
%     %plot(thetahist,col{i},'LineWidth',2)
%     hold on
%     grid on
%     xlim([0 90])
%     xlabel('Zenith angle [deg]', labelOpts{:})
%     ylabel('dN/d\theta [deg^{-1}]', labelOpts{:})
%
    phihist = t.PhiSTot;
    phitot = phitot + phihist;
    if i<7
        phitotew = phitotew+phihist;
    else
        phitotsn = phitotsn+phihist;
    end
    
%     figure(1)
%     subplot(2,1,1)
%     semilogy(phihist,col{i},'LineWidth',2)
%     %hold on
%     grid on
%     xlim([0 360])
%     ylabel('dN/d\phi [deg^{-1}]', labelOpts{:})
%     subplot(2,1,2)
%     %plot(phihist,col{i},'LineWidth',2)
%     hold on
%     grid on
%     xlim([0 360])
%     xlabel('Azimuth angle [deg]', labelOpts{:})
%     ylabel('dN/d\phi [deg^{-1}]', labelOpts{:})
    %
    xytot = xytot + t.XYTot;
    thetaphitot = thetaphitot + t.ThetaPhiTot;
    radtot = radtot + t.RadTot;
end
thetatot = thetatot./durtot;
phitot = phitot./durtot;

xrad = 10.^(0:0.1:9.9);
figure(45)
semilogx(xrad,radtot,'+k','MarkerSize',4,'LineWidth',1.5)
xlabel('Radius (m)')
grid on
sum(radtot)
sum(radtot(find(xrad>500)))
%pause

figure(1)
subplot(2,1,1)
semilogy(thetatot,'k','LineWidth',3)
hold on
semilogy(thetatotew,'k--','LineWidth',3)
grid on
xlim([0 90])
xlabel('Zenith angle (deg)', labelOpts{:})
ylabel('dN/d\theta (deg^{-1}day^{-1})', labelOpts{:})
subplot(2,1,2)
semilogy(phitot,'k','LineWidth',3)
hold on
semilogy(phitotew,'k--','LineWidth',3)
grid on
xlim([0 360])
xlabel('Azimuth angle (deg)', labelOpts{:})
ylabel('dN/d\phi (deg^{-1}day^{-1})', labelOpts{:})

figure(2)
subplot(2,1,1)
plot(thetatot,'k','LineWidth',2)
hold on
% plot(thetatotew,'b--','LineWidth',1)
% plot(thetatotsn,'r--','LineWidth',1)
grid on
xlim([0 90])
xlabel('Zenith angle (deg)', labelOpts{:})
ylabel('dN/d\theta (deg^{-1}day^{-1})', labelOpts{:})
subplot(2,1,2)
plot(phitot,'k','LineWidth',2)
hold on
% plot(phitotew,'b--','LineWidth',1)
% plot(phitotsn,'r--','LineWidth',1)
grid on
xlim([0 360])
xlabel('Azimuth angle (deg)', labelOpts{:})
ylabel('dN/d\phi (deg^{-1}day^{-1})', labelOpts{:})

% sum(thetatot(80:end))
% east=sum(phitot(70:90))
% west=sum(phitot(270:290))
% north=sum(phitot(350:360))+sum(phitot(1:10))
% (east+west+north)./sum(phitot)
%pause
disp(sprintf('All periods: %3.1e raw coincs, %3.1e after ConsCoinc; %3.1e after BadPulses',nrawtot,nconscoinctot,nbadtot))

figure(42)
polarplot3d(log10(thetaphitot+1));
hh = colorbar;
ylabel(hh, 'Event density rate (log_{10}(deg^{-2}day^{-1}))','Fontsize',12);

pause

figure(62)
phi=linspace(0,360,360);
theta=linspace(0,100,100);
%z = log10((thetaphitot+1)/durtot);
z = (thetaphitot+1)/durtot;
%h = surf(phi,theta,z);
Xq = linspace(min(phi),max(phi),length(phi)*10);
Yq = linspace(min(theta),max(theta),length(theta)*10);
[Xq,Yq] = meshgrid(Xq,Yq);
Zq = interp2(phi,theta,z,Xq,Yq,'linear');
h = surf(Xq,Yq,Zq);
xlabel('Azimuth (deg)')
ylabel('Zenith (deg)')
set(h,'LineStyle','None')
hh = colorbar;
ylabel(hh, 'Source density (log_{10}(deg^{-2}day^{-1}))','Fontsize',14);
ylim([0 90])
xlim([0 360])
view(0,-90)
%caxis([-2 2.5])
%return

figure(63)
% figure(39)
% subplot(2,2,1)
% errorbar(pstart:pstop,rawRate,eRawRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
% xlabel('Period ID', labelOpts{:})
% ylabel('Mean raw coinc rate [Hz]', labelOpts{:})
% hold on
% grid on
% % errorBar([pstart:pstop],selRate,eSelRate,'sr','MarkerSize',8,'MarkerFaceColor','r')
% subplot(2,2,3)
% errorbar(pstart:pstop,candRate,eCandRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
% xlabel('Period ID', labelOpts{:})
% ylabel('Candidate rate [day^{-1}]', labelOpts{:})
% hold on
% grid on
% subplot(1,2,2)
% plot(rawRate,candRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
% xlabel('Mean raw coinc rate [Hz]', labelOpts{:})
% ylabel('Candidate rate [day^{-1}]', labelOpts{:})
% hold on
% grid on

xytot=xytot+1;  % Add one event for better visibility
dens = log10(xytot'/440/durtot);
Xq = linspace(min(ew),max(ew),length(ew)*2);
%Xq = resample(ew,50,10);
Yq = linspace(min(sn),max(sn),length(sn)*2);
%Yq = resample(sn,50,10);
[Xq,Yq] = meshgrid(Xq,Yq);
Zq = interp2(ew,sn,dens,Xq,Yq,'linear');
figure(72)
h = surf(Xq,Yq,Zq);    
set(h,'LineStyle','None')
%surf(ew,sn,log10(xytot'/440),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
hold on
plot3(270,100,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(1416,-200,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(1477,-205,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(1765.9,426.4,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(1690.6,469.8,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(3347.6,356,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(3260.6,564.3,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(3414.3,560.2,max(Zq),'pm','MarkerFace', 'm', 'MarkerSize', 6 );
plot3(1490,-50,max(Zq),'sr','MarkerFace', 'r', 'MarkerSize', 8 );
xlabel('Easting (m)')
ylabel('Northing (m)')
zlabel('Source density (m$^{-2}$)')
h = colorbar;
ylabel(h, 'Source density (log_{10}(m^{-2}day^{-1}))','Fontsize',14);
xlim([-500,3500])
%xlim([min(ew),max(ew)])
ylim([-700 2000])
caxis([-5 0])
%ylim([min(sn),max(sn)])
view(0,90)
%axis equal

