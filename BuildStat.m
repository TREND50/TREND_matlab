function [ ] = BuildStat(runstart,runstop)
% Compute statistics for acquisition over a given PERIOD of time,
% based on log files & calib dsts.
% OMH 29/02/2012

SharedGlobals;

if ~exist('runstop')
    runstop = runstart;
end
looptime = 2^28*5e-9;
%trun = [2685:2890]; % March -> June 2011
%trun = [3000:3086]; % August -> October 2011
%trun = [3157:3256]; % November 2011 -> January 2012
%trun = [3337:3560]; % January -> February 2012
%trun = [3561:3710]; % February -> May 2012
trun = [runstart:runstop];
% Detectors IDs are hardcoded
Detectors = [101:140 148:158 175 176]; % Hardcoded
Scints = [139 175 176]; % Hardcoded
nDets = length(Detectors);
[c,id] = intersect(Detectors,Scints);
DetectorType = zeros(1,nDets);
DetectorType(id) = 1;
time_step_log = 2^28*5e-9; % Time to read a buffer; [s]
nDets = length(Detectors);
statdst = [MONITOR_PATH 'stat.mat'];

if fopen(statdst)<0
  tot = zeros(1,55);
  monit = zeros(1,54);
  up = zeros(1,54);
  good = zeros(1,54);
  crash = zeros(1,54);
  nspike = zeros(1,54);
  ntrig = zeros(1,54);
  t_t0 = zeros(1,54);
  t_t1 = zeros(1,54);
  t_loop = zeros(1,54);
  ngood = zeros(1,54);
  nsta = zeros(1,54);
  nbad = zeros(1,54);
  nsat = zeros(1,54);
  daqdur = zeros(1,2);
  k = 0;
else
  stat = load(statdst);
  tot = stat.tot;
  monit = stat.monit;
  up = stat.up;
  good = stat.good;
  crash = stat.crash;
  nspike = stat.nspike;
  ntrig = stat.ntrig;
  t_t0 = stat.t_t0;
  t_t1 = stat.t_t1;
  t_loop = stat.t_loop;
  ngood = stat.ngood;
  nsta = stat.nsta;
  nbad = stat.nbad;
  nsat = stat.nsat;
  %daqdur = stat.daqdur;
  k = size(tot,1);
end
nodatarun = [2576,2582,2585,3607,3721,3845,3862,3866]; % Runs for which no DAQ data was taken

for i=1:length(trun)
    if sum(ismember(nodatarun,trun(i)))>0
        disp(sprintf('No data in run %d, skipping it...',trun(i)))
        continue
    end
         
    %% Load log summary file
    filename = [MONITOR_PATH sprintf( 'summary_R%d.txt', trun(i) )]; 
    if fopen(filename)<0
        disp(sprintf('No run %d.',trun(i)))
        continue
    else
        if sum(ismember(tot(:,54),trun(i)))>0
            disp(sprintf('Run R%d already present in stat.mat, skipping it...',trun(i)))
            continue
        end
        disp(sprintf('Analysing run R%d...',trun(i)))
        fclose all;
    end

    summaryfile=load(filename);
    det = summaryfile(:,1);
    if trun(i)>=2685
        ns = summaryfile(:,2); % total number of spikes for each antenna
        nt = summaryfile(:,3); % total number of events recorded for each antenna
        if nt>ns
            disp 'Error in summary file on trigger & spike counts.'
            nt = 0;
            ns = 0;
        end
    else  % For Runs before 2685, counts in log files are not reliable.
        ns = 100*ones(1,length(det)); % total number of spikes for each antenna
        nt = 100*ones(1,length(det)); % total number of events recorded for each antenna        
    end
    tdeb = summaryfile(:,4); % time of 1st trigger for each antenna
    tend = summaryfile(:,5); % time of last trigger for each antenna
    if trun(i)>=3560 & size(summaryfile,2)==12 % new DAQ version
        tt0 = summaryfile(:,6); %[s]
        tt1 = summaryfile(:,7); %[s]
        tloop = summaryfile(:,8); %[s]
        ng = summaryfile(:,9); 
        nb = summaryfile(:,10);
        nst = summaryfile(:,11);
        nsa = summaryfile(:,12);
    else
        tt0 = zeros(1,length(det));
        tt1 = zeros(1,length(det));
        tloop = zeros(1,length(det));
        ng = zeros(1,length(det));
        nb = zeros(1,length(det));
        nst = zeros(1,length(det));
        nsa = zeros(1,length(det));
    end
    if size(find(tdeb>0),1) == 0
        disp 'Summary data corrupted. Skip this run.'
        continue
    end
    dur = (tend-tdeb)/60; % True acquisition time for each detector [minutes]
    
    if ismember(det,Scints)
        disp(sprintf('R%d (%3.1f mins) scints only... skiping it.',trun(i),dur(1)))
        fclose all;
        continue;
    end
    alert = find(dur<0);
    if length(alert)>0
        disp(sprintf('Error! Negative duration for %d antennas in R%d...\nCorrupted log files?. \nWill use max duration time (%3.1f minutes) for these antennas.',length(alert),trun(i),max(dur)))
        det(alert)
        pause
        dur(alert)=max(dur);
        continue
    end
    goodtime = find(abs(tdeb-max(tdeb))<1e3 & tdeb>0);  % Valid timing
    date = min(tdeb(goodtime));
    %
    [a idets idetm] = intersect(det,Detectors);
    k = k+1;  % k is index of valid (ie with antennas) run
    tot(k,idetm) = dur(idets);
    tot(k,55) = date;
    tot(k,54) = trun(i);
    monit(k,54) = trun(i);
    up(k,54) = trun(i);
    good(k,54) = trun(i);
    crash(k,54) = trun(i);
    nspike(k,54) = trun(i);
    ntrig(k,54) = trun(i);
    t_t0(k,54) = trun(i);
    t_t1(k,54) = trun(i);
    t_loop(k,54) = trun(i);    
    ngood(k,54) = trun(i);
    nsta(k,54) = trun(i);
    nsat(k,54) = trun(i);
    nbad(k,54) = trun(i);
    %
    %
    %% Load DAQ duration file (duration from log files)
    daqfile = [MONITOR_PATH sprintf( 'DAQduration_R%d.txt', trun(i) )]; 
    if fopen(daqfile)<0
        disp(sprintf('No DAQduration file for R%d.',trun(i)))
        daqdur(k,1) = -trun(i);
        daqdur(k,2) = max(dur)*60;  %s %Use duration from time files
    else
        DAQdurationfile = load(daqfile);
        unix_start = DAQdurationfile(1);
        unix_stop = DAQdurationfile(2);
        daqdur(k,1) = trun(i);
        daqdur(k,2) = unix_stop-unix_start;  %s
    end
    
    %% Load calib dst
    filename = [CAL_PATH sprintf('calib%d_1.mat',trun(i))];
    if fopen(filename)<0
        disp(sprintf('No calib dst available for R%d (ie no PSD data).',trun(i)))
        psd(k) = 0;
        %continue
    else
        fclose all;
        calib = load(filename);
        dets = calib.Struct.Setup.Antennas;
        nocal = calib.Struct.Calib.NoCalib;
        status = calib.Struct.Calib.Status;
        if isequal(dets,nocal)
            psd(k) = 0;
            disp(sprintf('Calib dst is empty for R%d (ie no PSD data).',trun(i)))
        else
            psd(k) = 1;
        end
    end
    % 
    %% Build summary matrixes
    for j = 1:nDets
        indsum = find(det==Detectors(j));  % index of detector in vectors from summary file.
        if size(indsum,1) == 0 % detector not in acquisition for this run -> skip
            disp(sprintf('Detector %d not in acquisition in this run.',Detectors(j))) 
            nspike(k,j) = 0;
            ntrig(k,j) = 0;
            monit(k,j) = 0;
            up(k,j) = 0;
            good(k,j) = 0;
            crash(k,j) = 0;
            t_t0(k,j) = 0;
            t_t1(k,j) = 0;
            t_loop(k,j) = 0;
            ngood(k,j) = 0;
            nsta(k,j) = 0;
            nbad(k,j) = 0;
            nsat(k,j) = 0;            
            continue
        end
        nspike(k,j) = ns(indsum);
        ntrig(k,j) = nt(indsum);
        t_t0(k,j) = tt0(indsum);
        t_t1(k,j) = tt1(indsum);
        t_loop(k,j) = tloop(indsum);
        nsta(k,j) = nst(indsum);     
        ngood(k,j) = ng(indsum);
        nbad(k,j) = nb(indsum);
        nsat(k,j) = nsa(indsum);     
        %
        if DetectorType(j) == 0  % Antennas: use PSD data to assess detector status.
            if psd(k) == 1
                indpsd = find(dets==Detectors(j));
                tmn = calib.Struct.Calib.Time(indpsd,:)/60;
                dtmn = diff(tmn);
                %ntot = length(tmn);  % not reliable to estimate timing...
                nmonit = length(find(dtmn>0));
                tstep = mean(dtmn(dtmn>0)); % step size of PSD measurements [mn]
                if isnan(tstep) 
                    continue
                end
                nup = 0;
                ngoo = 0;
                ntot = ceil(tot(k,j)/tstep);
                if nmonit > ntot+5
                    disp(sprintf('Antenna %d: Monitored time (%3.1f mins) exceeds total run time (%3.1f mins)!',Detectors(j),nmonit*tstep,tot(k,j)))
                    crash(k,j) = (nmonit-ntot)*tstep;  % If PSD duration is larger than acquisition, this means that DAQ crashed.
                    nmonit = ntot;
                end
                evo = status(indpsd,1:nmonit);    
                nneg = sum(evo==-1);
                if nneg==length(evo)  % No status for this antenna
                    nmonit = 0;
                elseif nneg>0  % Some status values are negative
                    disp 'Error! OBSOLETE format for calib dst (status = -1)'
                    pause
                else % No negative values
                    for l = 1:nmonit
                        st = dec2bin(evo(l));
                        stv = int16(sscanf(st,'%1d'))'; %LSB last
                        stv=fliplr(stv); %LSB first
                        if stv(1)==0
                          nup = nup + 1;  % Antenna not dead
                          if length(stv)==1 | stv(2)==0  % Antenna not bumpy
                              ngoo = ngoo+1;
                          end
                        end
                    end
                end        
                monit(k,j) = nmonit*tstep; % time for which status evaluation is possible (psd running) [min]
                up(k,j) = nup*tstep; % fft is here
                good(k,j) = ngoo*tstep;  % fft is good  
            end
        else  % Scints: no PSD data, have to rely on spike rate to assess detector state
            filename = sprintf( 'R%06d_A%04d_log.txt', trun(i), Detectors(j) ); 
            filename = [LOG_PATH, filename];
            if fopen(filename)>0  % This scint was in the DAQ
                disp(sprintf('Analysing file %s',filename))
                logfile=load(filename);
                ts = logfile(:,1); % Time at begining of reading loop
                irq =  logfile(:,2); % Buffer index
                strig = find(ts>0);
                %tot(k,j) = length(ts)*time_step_log*2/60*(dur(indsum)>0); %[min]
                tot(k,j) = irq(end)*time_step_log/60; %[min]
                monit(k,j) = tot(k,j); %[min]
                up(k,j) = 0;
                good(k,j) = 0;
                if size(strig,1)==0
                    disp(sprintf('No data recorded on detector %d.',Detectors(j)));
                else  % data recorded
                    t_firsttrig = strig(1)*time_step_log*2;
                    %parity = [diff(irq)' 0]';
                    if ts(1)>0
                        tmn = (ts-ts(1))/60;  % Time in minutes since run start... NOT GOOD IF NO TRIG on 1st line...
                    else
                        %tmn = t/60/2;
                        %t = t_firsttrig + irq.*time_step_log*2+parity.*time_step_log;
                        %tmn = t/60;
                        tmn = t_firsttrig + irq(strig(1)).*time_step_log/60;
                    end
                    nns = logfile(:,3); % number of spikes in this half buffer
                    nnt = logfile(:,4); % number of events recorded on disc in this half buffer 
                    rate = logfile(:,5); % Spike rate [Hz]                    
                    if ~exist('tstep')
                        tstep = 10;
                    end
                    disp(sprintf('Acquisition duration = %3.1f mn. Using time step = %3.1f mn',tmn(end),tstep))
                    for l=1:round(tmn(end)/tstep)
                      sel = find(tmn>=(l-1)*tstep & tmn<l*tstep);  % time slice
                      if sum(nns(sel)) > 0 
                          up(k,j) = up(k,j)+tstep;
                          rhz = sum(nns(sel))/tstep/60;
                          if rhz>5 & rhz<500  % Scintillator in "good" state if spike rate between 10 & 200Hz.
                              good(k,j) = good(k,j)+tstep;
                          end
                      end
                    end
                end
            else
                disp(sprintf('No file %s, skipping it.',filename))
            end          
        end % use log files
    end
end

tot = sortrows(tot,54);
monit = sortrows(monit,54);
up = sortrows(up,54);
good = sortrows(good,54);
crash = sortrows(crash,54);
ntrig = sortrows(ntrig,54);
nspike = sortrows(nspike,54);
t_t0 = sortrows(t_t0,54);
t_t1 = sortrows(t_t1,54);
t_loop = sortrows(t_loop,54);
nsta = sortrows(nsta,54);
ngood = sortrows(ngood,54);
nbad = sortrows(nbad,54);
[nsat ind] = sortrows(nsat,54);
%daqdur = daqdur(ind,:);  % Cannot use sortrows(daqdur,1) because run without DAQDuration files have runnumbers set to -runnumber
%
[a indrun] = intersect(tot(:,54),trun);  % Compute list of indexes in matrix for runs selected
for i = 1:length(indrun)
    irun = indrun(i);
    selbad = find(ntrig(irun,:)./nspike(irun,:)>1);
    for j=1:length(selbad)
        jbad = selbad(j);
        disp(sprintf('Warning! R%d Antenna %d: Wrong ratio NT1/NT0 = %d/%d = %3.2f pc',tot(irun,54),Detectors(jbad),ntrig(irun,jbad),nspike(irun,jbad),ntrig(irun,jbad)/nspike(irun,jbad)*100))
        %nspike(irun,jbad) = 0;
        %CleanStat(tot(irun,54),tot(irun,54))
        pause
    end
end

disp '************ No save at present ***************'
%save(statdst, 'tot', 'monit', 'up', 'good', 'crash', 'ntrig', 'nspike','t_t0','t_t1','t_loop','ngood','nbad','nsat','nsta','daqdur') 

%% Stats for this set
%[yi, mi, di, hi, mni, si] = UnixSecs2Date(tot(:,55));
%dc = datenum(yi,mi,di,hi,mni,si);
for j = 1:nDets
    albl{ j } = num2str( Detectors(j) );
    dtot(j) = sum(tot(indrun,j))/60/24;
    dmonit(j) = sum(monit(indrun,j))/60/24;
    dup(j) = sum(up(indrun,j))/60/24;
    dgood(j) = sum(good(indrun,j))/60/24;
    dcrash(j) = sum(crash(indrun,j))/60/24;
    nt0(j) = sum(nspike(indrun,j));
    nt1(j) = sum(ntrig(indrun,j));
    ssta(j) = sum(nsta(indrun,j));
    sgood(j) = sum(ngood(indrun,j));
    sbad(j) = sum(nbad(indrun,j));
    ssat(j) = sum(nsat(indrun,j));
    in = find(t_t0(:,j)>0 & t_t0(:,j)<looptime);
    in = intersect(indrun,in);
    avt0(j) = mean(t_t0(in,j))*100;
    avt1(j) = mean(t_t1(in,j))*100;    
    avloop(j) = mean(t_loop(in,j))*100;
end
% Do not consider treatùent time 
avt0(DetectorType==1) = 0;
avt1(DetectorType==1) = 0;
avloop(DetectorType==1) = 0;

%DAQRunningDays = sum(daqdur(indrun,2))/3600/24;  % Acquisition running (log files written) [days]
DAQRunningDays = 0
totloop = sgood+sbad+ssta;
rmonit = dmonit./dtot*100;
rup = dup./dmonit*100;
rgood = dgood./dtot*100;
%rcrash = dcrash./dtot*100;
rt1 = nt1./nt0*100;
ddeb = tot(indrun(1),55);
durlast = max(tot(indrun(end),1:53))*60; % duration of last run [secs]
dend = tot(indrun(end),55) + durlast; % date of end of last run [unix secs]
durd = (dend-ddeb)/3600/24; 
rdump = sbad./totloop*100;  % Duumped buffers
rsat = ssat./totloop*100; % Satturated buffers
rsta = ssta./totloop*100; % Stalled buffers
rbad = (ssta+sbad)./totloop; % Stalled+Dumped
rdump(isnan(rdump)) = 0;
rsat(isnan(rsat)) = 0;
rsta(isnan(rsta)) = 0;
rbad(isnan(rbad)) = 0;
r_t0 = avt0/looptime;
r_t1 = avt1/looptime;
r_loop = avloop/looptime;
[yd md dd hd mnd]=UnixSecs2Date(ddeb);
[ye me de he mne]=UnixSecs2Date(dend);
disp(sprintf('\n*** Runs %d-%d (%02d/%02d/%d - %02d/%02d/%d : %3.1f true days - %3.1f DAQ days)\n',tot(indrun(1),54),tot(indrun(end),54),dd,md,yd,de,me,ye,durd,DAQRunningDays))
disp '   Id     Total[days] Monitored[%] Up[%]   Good[%]   NbT1[G]   T1/T0[%]'
[Detectors' dtot' rmonit'  rup' rgood' nt1'/1e9 rt1']

%% Plot
figure(3)
tit = sprintf('Stats R%d-R%d (%02d/%02d/%d - %02d/%02d/%d : %3.1f true days - %3.1f DAQ days)',tot(indrun(1),54),tot(indrun(end),54),dd,md,yd,de,me,ye,durd,DAQRunningDays);
set(3,'Name',tit,'NumberTitle','off')
subplot(4,1,1)
plot(1:nDets,dtot,'ks','MarkerFaceColor','w','MarkerSize',9)
grid on
hold on
plot(1:nDets,dmonit,'ks','MarkerFaceColor','k','MarkerSize',9)
plot(1:nDets,dup,'gs','MarkerFaceColor','g','MarkerSize',8)
l = line([0 nDets+1],[DAQRunningDays DAQRunningDays]);
set(l,'Color','k','LineWidth',2)
xlim([0 nDets+2])
grid on;
set( gca, 'FontSize', 10, 'FontWeight', 'bold' );
set( gca, 'XTick', 1:nDets, 'XTickLabel', albl );
ylabel( 'Live time [days]' );
xlabel( 'Detector Id' );
legend('time in DAQ','time with PSD','time up (PSD)')
title(tit)
%
subplot(4,1,2)
semilogy(1:nDets,nt0+1,'ks','MarkerFaceColor','k','MarkerSize',8)
hold on
semilogy(1:nDets,nt1+1,'gs','MarkerFaceColor','g','MarkerSize',8)
xlim([0 nDets+2])
grid on;
set( gca, 'FontSize', 10, 'FontWeight', 'bold' );
set( gca, 'XTick', 1:nDets, 'XTickLabel', albl );
ylabel( 'Nevents' );
xlabel( 'Detector Id' );
legend('N T0','N T1')
%
subplot(4,1,3)
plot(1:nDets,r_t0,'bs','MarkerFaceColor','b','MarkerSize',8)
hold on
grid on
plot(1:nDets,r_t1,'gs','MarkerFaceColor','g','MarkerSize',8)
plot(1:nDets,r_loop,'ks','MarkerFaceColor','k','MarkerSize',8)
xlim([0 nDets+2])
ylim([0 110])
set( gca, 'FontSize', 10, 'FontWeight', 'bold' );
set( gca, 'XTick', 1:nDets, 'XTickLabel', albl );
ylabel( 'Anal. time [%]' );
xlabel( 'Detector Id' );
legend('mean T0 time', 'mean T1 time','mean total')
%
subplot(4,1,4)
plot(1:nDets,rsat,'ks','MarkerFaceColor','k','MarkerSize',8)
hold on
grid on
plot(1:nDets,rdump,'rs','MarkerFaceColor','r','MarkerSize',8)
xlim([0 nDets+2])
ylim([0 110])
plot(1:nDets,rsta,'rs','MarkerFaceColor','w','MarkerSize',8)
xlim([0 nDets+2])
set( gca, 'FontSize', 10, 'FontWeight', 'bold' );
set( gca, 'XTick', 1:nDets, 'XTickLabel', albl );
ylabel( 'Bad & sat buffers [%]' );
xlabel( 'Detector Id' );
legend('% max trig rate','% dumped buffers','% stalled DAQ')

figure(17)
[yi, mi, di, hi, mni, si] = UnixSecs2Date(tot(indrun,55));
d = datenum(yi,mi,di,hi,mni,si);
plot(d,tot(indrun,54),'+')
datetick('x','mm/yyyy')

DAQRunningDays
goodAnts = (DetectorType==0) & isfinite(dtot); 
DAQTotDays = max(dtot(goodAnts))
BadFrac = mean(rbad(goodAnts))
DAQliveDays = DAQTotDays*(1-BadFrac)


