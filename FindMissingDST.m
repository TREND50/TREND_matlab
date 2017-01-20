function [] = FindMissingDST(runstart,runstop)
% Find runs present in stat.mat but with no "Candidate" DST.
% OMH 04/02/2013

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
j = 0;
missing = [0 0];
present = [0 0];
for i = runstart:runstop

    %% Get BuildStat results
    irun = find(ttot(:,54)==i);
    if size(irun,1) == 0
        disp(sprintf('Run %d not in %s. Skip.',i,statdst))
        continue
    end
    t_tot = ttot(irun,DetectorType==0)/60; % DAQ physical duration [h]
    
    j = j+1;
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
    disp(sprintf('R%d: %3.1f h',i,max(t_tot)))
    if nbiter==0
        display(sprintf('No dst found for run %d. Resseting live time computation.',i))
        dTot(j) = 0;
        dLive(j) = 0;
        EffSurf(j) = 0;
        CheckTot(j) = 0;
        missing(end+1,1) = i;
        missing(end,2) = max(t_tot);
        continue 
    else
        present(end+1,1) = i;
        present(end,2) = max(t_tot);
    end;
    %%
end
missing(1,:) = [];
filename = sprintf('missing_R%dR%d.mat',runstart,runstop);
save(filename,'missing')

inds = find(missing(:,2)<1);
indl = find(missing(:,2)>=1);
ds = sum(missing(inds,2));
dl = sum(missing(indl,2));
dp = sum(present(:,2));

disp(sprintf('R%d - R%d:',runstart,runstop))
disp(sprintf('%d runs present. Total = %3.1f h',length(present),dp))
disp(sprintf('%d missing runs with duration less than 1h. Total = %3.1f h',length(inds),ds))
disp(sprintf('%d missing runs with duration >= 1h. Total = %3.1f h:',length(indl),dl))
[missing(indl,1) missing(indl,2)]

