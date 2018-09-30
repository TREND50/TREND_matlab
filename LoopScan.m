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
timetot = [];
runtot = [];
if EW==1
    pstart = 1;
    pstop = 8;
else
    pstart = 9;
    pstop = 10;
end

down
for i=pstart:pstop
   [thp php mult rad ramp run a b t] = ScanCandidateDST(i,down);
   sel = find(php>0 & thp>0 & thp<=80);
   %sel = find(thp<=80 & php>=270 & php<=360);
   run = run(sel);
   t = t(sel);
   thp = thp(sel);
   php = php(sel);
   mult = mult(sel);
   rad = rad(sel);
   ramp = ramp(sel);
   runtot = [runtot run];
   timetot = [timetot t];
   thtot = [thtot thp];
   phtot = [phtot php];
   multot = [multot mult];
   radtot = [radtot rad];
   ramptot = [ramptot ramp];
   
   disp(sprintf('%d candidates selected in period %d.',length(thp),i))
   %pause
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
errorbar(tht,N,sqrt(N), 'sk-','MarkerFaceColor','k','LineWidth',2 );
%errorbar(tht,N/N(7)*0.3425,sqrt(N)/max(N), 'sb-','MarkerFaceColor','b','LineWidth',2 );
grid on
xlim([0 90])
%ylim([0 1.9])
xlabel('Zenith angle [deg]', labelOpts{:})
ylabel('dN/d\theta [10deg^{-1}]', labelOpts{:})
disp ' '

figure(13)
hold on
N = [];
for i=1:18
    minedge = (i-1)*20;
    maxedge = i*20;
    minsel = mod(minedge-90,360);
    maxsel = mod(maxedge-90,360);
    if minsel>maxsel
       sel = find( (phtot>=minsel | phtot<maxsel) & thtot<80);
    else
        sel = find(phtot>=minsel & phtot<maxsel & thtot<80);
    end
    %pht(i) = mean(phtot(sel));
    pht(i) = (maxedge+minedge)/2;
    N(i) = length(sel);
end    
errorbar(pht,N,sqrt(N), 'sr-','MarkerFaceColor','r','LineWidth',2 );
xlabel('Azimuth angle [deg]', labelOpts{:})
ylabel('dN/d\phi [40deg^{-1}]', labelOpts{:})
grid on
xlim([0 360])

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
prepareSkyplot(4);
h = polar( phtot*DEG2RAD(1), thtot, 'k*');
set(h, 'MarkerFaceColor', 'g' );  
if EW==0 % Add simulated data
  simNS = load("NSsim.txt");
  thetaSim = simNS(:,1);
  phiSim = mod(simNS(:,2)-90,360);
end
polar( phiSim*DEG2RAD(1), thetaSim, 'r*');
t = length(phiSim);
n = length(find(phiSim<90 | phiSim>270));
w = length(find(phiSim<180));
 
disp ' '
disp(sprintf('Total nb of sim candidates = %d.',t))
disp(sprintf('North/South = %d/%d (%3.1f pc excess)',n,t-n,n/t*100))
disp(sprintf('West/East = %d/%d (%3.1f pc excess)',w,t-w,w/t*100))

SmoothSkyplot(thtot,phtot);

t = length(phtot);
n = length(find(phtot<90 | phtot>270));
w = length(find(phtot<180));

nw = length(find(phtot>0 & phtot<90));
sw = length(find(phtot>90 & phtot<180));
se = length(find(phtot>180 & phtot<270));
ne = length(find(phtot>270 & phtot<360));
 
disp(sprintf('Total nb of candidates = %d.',t))
disp(sprintf('North/South = %d/%d (%3.1f pc excess)',n,t-n,n/t*100))
disp(sprintf('West/East = %d/%d (%3.1f pc excess)',w,t-w,w/t*100))
s = w-(t-w)
s/t*100

fclose all;
filename = 'candidates.txt';
fid = fopen( filename, 'w' );
for i=1:length(timetot)
    fprintf(fid,'%20d ',timetot(i));
    fprintf(fid,'%3d ',runtot(i));
    fprintf(fid,'%3.2f ',thtot(i));
    fprintf(fid,'%3.2f ',phtot(i));
    fprintf( fid, '\n' );
end
fclose(fid);

