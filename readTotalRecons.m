function [] = readTotalRecons()
% OMH 05/06/2017


SharedGlobals;
pstart = 1;
pstop = 8;
dayToSecs = 3600*24;
dur = [22.8 42.2 44.1 49.6 23.2 79.2 26.2 28.6];   % Taken from python monitoring tool for period >=6 (/sps/hep/TREND/readSummaryProd)
sum(dur)
cands = [44 94 87 44 54 205 24 22];
candRate = cands(pstart:pstop)./dur(pstart:pstop);
eCandRate = sqrt(cands(pstart:pstop))./dur(pstart:pstop);


col = {'k','b','g','r','m','c','y','b--',};
phitot = zeros(1,360);
thetatot = zeros(1,100);
nrawtot = 0;
nconscoinctot = 0;
nbadtot = 0;
rawRate = zeros(1,pstop-pstart+1);
selRate = zeros(1,pstop-pstart+1);
eRawRate = zeros(1,pstop-pstart+1);
eSelRate = zeros(1,pstop-pstart+1);


for i=pstart:pstop
    totdst = sprintf('TotalRecons_Period%d_102014.mat',i);
    totdst = [CAND_PATH totdst];
    t = load(totdst);
    runs = t.Runs;
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
    thetahist = t.ThetaPTot;
    thetatot = thetatot + thetahist;
    figure(2)
    subplot(2,1,1)
    semilogy(thetahist,col{i},'LineWidth',2)
    %hold on
    grid on
    xlim([0 90])
    ylabel('dN/d\theta [deg^{-1}]', labelOpts{:})
    subplot(2,1,2)
    %plot(thetahist,col{i},'LineWidth',2)
    hold on
    grid on
    xlim([0 90])
    xlabel('Zenith angle [deg]', labelOpts{:})
    ylabel('dN/d\theta [deg^{-1}]', labelOpts{:})
    %
    phihist = t.PhiPTot;
    phitot = phitot + phihist;
    figure(1)
    subplot(2,1,1)
    semilogy(phihist,col{i},'LineWidth',2)
    %hold on
    grid on
    xlim([0 360])
    ylabel('dN/d\phi [deg^{-1}]', labelOpts{:})
    subplot(2,1,2)
    %plot(phihist,col{i},'LineWidth',2)
    hold on
    grid on
    xlim([0 360])
    xlabel('Azimuth angle [deg]', labelOpts{:})
    ylabel('dN/d\phi [deg^{-1}]', labelOpts{:})
end

figure(1)
subplot(2,1,1)
semilogy(phitot,'k','LineWidth',3)

subplot(2,1,1)
semilogy(phitot,'k','LineWidth',3)
grid on
xlim([0 360])
ylabel('dN/d\phi [deg^{-1}]', labelOpts{:})
subplot(2,1,2)
plot(phitot,'k','LineWidth',3)
figure(2)
subplot(2,1,1)
semilogy(thetatot,'k','LineWidth',3)
grid on
xlim([0 90])
ylabel('dN/d\theta [deg^{-1}]', labelOpts{:})
subplot(2,1,2)
plot(thetatot,'k','LineWidth',3)

disp(sprintf('All periods: %3.1e raw coincs, %3.1e after ConsCoinc; %3.1e after BadPulses',nrawtot,nconscoinctot,nbadtot))

figure(39)
subplot(2,2,1)
errorBar([pstart:pstop],rawRate,eRawRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
xlabel('Period ID', labelOpts{:})
ylabel('Mean raw coinc rate [Hz]', labelOpts{:})
hold on
grid on
% errorBar([pstart:pstop],selRate,eSelRate,'sr','MarkerSize',8,'MarkerFaceColor','r')
subplot(2,2,3)
errorBar([pstart:pstop],candRate,eCandRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
xlabel('Period ID', labelOpts{:})
ylabel('Candidate rate [day^{-1}]', labelOpts{:})
hold on
grid on
subplot(1,2,2)
plot(rawRate,candRate,'sk','MarkerSize',8,'MarkerFaceColor','k')
xlabel('Mean raw coinc rate [Hz]', labelOpts{:})
ylabel('Candidate rate [day^{-1}]', labelOpts{:})
hold on
grid on
