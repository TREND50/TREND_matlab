function [] = LoopScan(EW,down)
% Loop ScanCandidateDST on all periods.
% Produce smoothed skyplot
% OMH 27/06/2013

SharedGlobals;
if ~exist('down')  
  down = 0;
end
thtot = [];
phtot = [];
multot = [];
radtot = [];
ramptot = [];

if EW==1
    pstart = 1;
    pstop = 8;
else
    pstart = 9;
    pstop = 10;
end

down
for i=pstart:pstop
   [thp php mult rad ramp] = ScanCandidateDST(i,down);
   sel = find(php>0 & thp>0 & thp<=80);
   %sel = find(thp<=80 & php>=270 & php<=360);
   thp = thp(sel);
   php = php(sel);
   mult = mult(sel);
   rad = rad(sel);
   ramp = ramp(sel);
   thtot = [thtot thp];
   phtot = [phtot php];
   multot = [multot mult];
   radtot = [radtot rad];
   ramptot = [ramptot ramp];
   
   disp(sprintf('%d candidates selected in period %d.',length(thp),i))
   pause
end

figure(5)
subplot(2,1,1)
hist(thtot,15)
xlim([0 90])
xlabel('Zenith angle [deg]',labelOpts{:})
subplot(2,1,2)
hist(phtot,60)
xlim([0 360])
xlabel('Azimuth angle [deg]',labelOpts{:})

figure(12)
hold on
for i=1:9
    minedge = (i-1)*10;
    maxedge = i*10;
    sel = find(thtot>=minedge & thtot<maxedge & thtot<80 & thtot>0);
    tht(i) = mean(thtot(sel));
    N(i) = length(sel);
end
errorbar(tht,N/max(N),sqrt(N)/max(N), 'sb-','MarkerFaceColor','b' );
%errorbar(tht,N/N(7)*0.3425,sqrt(N)/max(N), 'sb-','MarkerFaceColor','b','LineWidth',2 );
grid on
xlim([0 90])
%ylim([0 1.9])
xlabel('Zenith angle [deg]', labelOpts{:})
ylabel('Fraction of showers detected', labelOpts{:})
disp ' '

figure(13)
hold on
N = [];
for i=1:36
    minedge = (i-1)*10+0;
    maxedge = i*10+0;
    sel = find(phtot>=minedge & phtot<maxedge & thtot<80);
    pht(i) = mean(phtot(sel));
    N(i) = length(sel);
end    
errorbar(pht,N/max(N),sqrt(N)/max(N), 'sk-','MarkerFaceColor','k' );
grid on
xlim([0 361])

figure(8)
hist(multot(multot>=5))
%xlim([5 max(multot)])
xlabel('Antenna multiplicity',labelOpts{:})


figure(200)
hist(ramptot,300)
xlabel('Amplitude Ratio', labelOpts{:})
grid on

figure(300)
subplot(2,1,1)
hist(log10(radtot),100)
xlabel('log10(Radius) [m]', labelOpts{:})
subplot(2,1,2)
semilogy(thtot,radtot,'k+')
grid on
xlabel('\theta [deg]', labelOpts{:})
ylabel('log10(Radius) [m]', labelOpts{:})

%pause
figure(4)
PrepareSkyPlot(4);
h = polar( phtot*DEG2RAD(1), thtot, 'k*');
set(h, 'MarkerFaceColor', 'g' );  
SmoothSkyplot(thtot,phtot);

t = length(phtot);
n = length(find(phtot<90 | phtot>270));
w = length(find(phtot<180));



nw = length(find(phtot>0 & phtot<90))
sw = length(find(phtot>90 & phtot<180))
se = length(find(phtot>180 & phtot<270))
ne = length(find(phtot>270 & phtot<360))
 

disp ' '
disp(sprintf('Total nb of candidates = %d.',t))
disp(sprintf('North/South = %d/%d (%3.1f pc excess)',n,t-n,n/t*100))
disp(sprintf('West/East = %d/%d (%3.1f pc excess)',w,t-w,w/t*100))
s = w-(t-w)
s/t*100
[thtot' phtot']