function [] = ComputeAcceptanceTotal(runstart,runstop)
% Cumulated acceptance in h.m2
% OMH 30/11/2012

SharedGlobals;
time_step_log = 2^28*5e-9;

Detectors = [101:140 148:158 175 176]; % Hardcoded
Scints = [139 175 176]; % Hardcoded
nDets = length(Detectors);
[c,id] = intersect(Detectors,Scints);
DetectorType = zeros(1,nDets);
DetectorType(id) = 1;  % Scints: id=1; Ants: id=0
Ants = Detectors(DetectorType==0);

statdst = [MONITOR_PATH 'stat.mat'];
stat = load(statdst);
ttot = stat.tot;

Runs = [];
dLive = [];
EffSurf = [];
CheckTot = [];
sensr = [];
nAnts = [];
j = 0;
for i = runstart:runstop

    %% Get BuildStat results
    irun = find(ttot(:,54)==i);
    if size(irun,1) == 0
        disp(sprintf('Run %d not in %s. Skip.',i,statdst))
        continue
    end
    j = j+1;
    
    % Live time
    t_tot = stat.tot(irun,DetectorType==0)/60/24; % DAQ physical duration [days]
    if i>3560  % New format for log files
        ngood = stat.ngood(irun,DetectorType==0);
        nbad = stat.nbad(irun,DetectorType==0)+stat.nsta(irun,DetectorType==0);
        nloops = ngood+nbad;
        check = nloops*time_step_log/3600/24-t_tot;
        rgood = ngood./nloops;  % Fraction of "good" buffers
    else  
        rgood = 0.8*ones(1,length(Ants));  % Assuming a 80% data quality from R3562-R3659 set
        check = zeros(1,length(Ants));
    end
    tlive = rgood.*t_tot; % Live DAQ time per antenna
    tlive(isnan(tlive))=0;
    Runs(j) = i;
    dTot(j) = max(t_tot);
    dLive(j) = max(tlive);
    if dLive(j)>0
        effAnt = tlive/dLive(j);  % Live DAQ time expressed as a ratio of max live time.
    end
    % Valid antennas
    nt0 = stat.nspike(irun,DetectorType==0);
    nt1 = stat.ntrig(irun,DetectorType==0);
    iout = find(nt1<100);  % Not enough triggers: antenna can be considered as dead
    iin = find(nt1>=100);
    nAnts(j) = length(iin);
%     rt1 = log10(nt1(iin)./nt0(iin));
%     pt1 = (log10(nt1(iin))-log10(mean(nt1(iin))))./log10(mean(nt1(iin)))
%     Ants(iin(pt1<-0.5))
%     Ants(iout)
%     figure; 
%     hist(rt1)
%     figure; 
%     hist(pt1,50)
    effAnt(iout==1) = 0; 
    %
    %EffSurf(j) = sum(effAnt)/length(Ants);
    if nAnts(j) == 0
        EffSurf(j) = 0;    
    else
        EffSurf(j) = sum(effAnt(iin))/nAnts(j);
    end
    CheckTot(j) = max(abs(check));
    disp(sprintf('Run %d: %3.2f live days, effective area ratio = %3.2f pc',i,dLive(j),EffSurf(j)*100))
    
    %continue
    
    %% Get number of sub dsts
    stopflag=0;
    nbiter=1;
    while stopflag==0
        filename = [DST_PATH sprintf('dst%d_%d.mat',i,nbiter)];
        fd=fopen(filename);
        if fd~=-1
            nbiter=nbiter+1;
            fclose(fd);
        else
            stopflag=1;
        end;
    end;
    nbiter=nbiter-1;
    
    if nbiter==0
        display(sprintf('No dst found for run %d. Resseting live time computation.',i))
        dTot(j) = 0;
        dLive(j) = 0;
        EffSurf(j) = 0;
        CheckTot(j) = 0;
        sensr(j,:) = zeros(1,360);
        continue
    end;
    
    %% Loop on sub dsts
    for k=1:nbiter        
%         filename = [DST_PATH sprintf('dst%d_%d.mat',i,k)];
%         display(sprintf('Loading run %d_%d',i,k))
%         clear d
%         d = load(filename);
        %
        filename = [ACC_PATH sprintf('Acceptance_R%d_%d.mat',i,k)];
        fd=fopen(filename);
        if fd~=-1 % Acceptance was already computed
            a = load(filename);
            sens(k,:) = a.sensall;
            dur(k) = a.dur;
            clear a;
        else
          disp 'Now calling ComputeAcceptance...'
          [coco sens(k,:) dur(k) ]= ComputeAcceptance(i,k);
        end
    end
    sens(~isfinite(sens)) = 0;
    sens(isnan(sens)) = 0;
    fclose all;
    sensr(j,:) = zeros(1,360);
    for k = 1:nbiter
        sensr(j,:) = sensr(j,:)+sens(k,:)*dur(k);  % Compute run acceptance by weighting subdsts with duration  
    end
    sensr(j,:) = sensr(j,:)./sum(dur);
end
sensr(~isfinite(sensr)) = 0;
sensr(isnan(sensr)) = 0;

%%
senstot = zeros(1,360);
for j = 1:length(dLive)
    senstot = senstot+sensr(j,:)*dLive(j)*EffSurf(j)*1.5*nAnts(j)/50;
end

%% Save to file
filename = sprintf('TotalAcceptance_R%dR%d.mat',runstart,runstop);
save([ACC_PATH filename],'Runs','nAnts','dTot','dLive','EffSurf','CheckTot','senstot','sensr');

% % Plot
% figure(11);
% set(11,'Name','Total acceptance','NumberTitle','off')
% plot(1:360,senstot,'k','LineWidth',2)
% xlabel('Azimuth angle [deg]',labelOpts{:})
% ylabel('Total acceptance [days x km2]',labelOpts{:})
% xlim([0 360])
% grid on
% hold on
%
for i=1:length(Runs)
    disp(sprintf('R%d %3.2f',Runs(i),dLive(i)))
end
sum(dTot)
sum(dLive)
sum(EffSurf.*dLive)/sum(dLive)
sum(CheckTot)
