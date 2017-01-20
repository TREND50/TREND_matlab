function [] = CandidateAnalysis_NS(nrun,dstf)
% Select candidates in NS polar data
% Writes it to dst.
% OMH 11/06/2013

SharedGlobals;
cutsettings.ThetaCut = 85;
cutsettings.RCut = 1000;  % % Mini radius [m]
cutsettings.BaryCut = 500; % Mean distance to coinc barycenter [m]
cutsettings.TimeCut = 30; % Time window [seconds]
cutsettings.DirTimeCut = 10; % Time window for same direction [minutes]
cutsettings.MaxAnt = 10;
cutsettings.DirMaxAnt = 3;
cutsettings.PhiCut = 10;  % Same direction [deg]
cutsettings.AntRatioCut = 0.5; % Max ratio of antennas in common with neighbouring events 
cutsettings.AmpRatioCut = 1.5;


%% Load valid Candidates DST
periods = zeros(10,2);
periods(6,:) = [3562 3733];
periods(9,:) = [4444 4844];  %NS polar
periods(10,:) = [4845 5070]; %NS polar + DAQ upgrade
periodID = find(periods(:,1)<=nrun & periods(:,2)>=nrun);
if size(periodID,1)==0
    disp(sprintf('Run %d does not belong to an identified period! Abort.',nrun))
    return
end
filename = sprintf('Candidates_Period%d.mat',periodID);
disp(sprintf('Now loading DST %s...',filename))
if fopen(filename)>0
    % load dst
    c = open(filename);
    CandidateRun = c.CandidateRun;
    CandidateCoinc = c.CandidateCoinc;
    CandidateTime = c.CandidateTime;
    CandidateThetaP = c.CandidateThetaP;
    CandidatePhiP = c.CandidatePhiP;
    CandidateThetaS = c.CandidateThetaS;
    CandidatePhiS = c.CandidatePhiS;
    CandidateRadius = c.CandidateRadius;
    CandidateAntennas = c.CandidateAntennas;
    CandidateNeighbourgs = c.CandidateNeighbourgs;
    CandidateDirNeighbourgs = c.CandidateDirNeighbourgs;
    %
    %Clear data for this run
    thisRun = find(CandidateRun==nrun);
    CandidateThetaS(thisRun) = [];
    CandidatePhiS(thisRun) = [];
    CandidateRadius(thisRun) = [];
    CandidateAntennas(thisRun) = [];
    CandidateDirNeighbourgs(thisRun) = [];
    CandidateNeighbourgs(thisRun) = [];
    CandidatePhiP(thisRun) = [];
    CandidateThetaP(thisRun) = [];
    CandidateTime(thisRun) = [];
    CandidateCoinc(thisRun) = [];
    CandidateRun(thisRun) = [];
else
    CandidateRun=[];
    CandidateCoinc=[];
    CandidateTime=[];
    CandidateThetaP=[];
    CandidatePhiP=[];
    CandidateThetaS=[];
    CandidatePhiS=[];
    CandidateRadius=[];
    CandidateTotal=0;
    CandidateAntennas={};
    CandidateNeighbourgs={};
    CandidateDirNeighbourgs={};
end
disp 'Done.'

%% Get number of sub dsts
stopflag=0;
nbiter=1;
while stopflag==0
    filename = [DST_PATH sprintf('dst%d_%d_light.mat',nrun,nbiter)];
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
    display(sprintf('No dst found for run %d',nrun))
else
    display(sprintf('%d dst(s) found for run %d.',nbiter,nrun))
end;

%% Loop on sub dsts
meta = 1;
while meta<=nbiter

    %% Load dst
    if ~exist('dstf')
        %dstname = [DST_PATH sprintf(dst_filename,nrun,meta)];
        dstname = [DST_PATH sprintf('dst%d_%d_light.mat',nrun,meta)];
        disp(sprintf('Loading dst %d for run %d...',meta,nrun))
        dst = load(dstname);
        disp 'Done.'
    else
        disp 'Using dst passed as argument.'
        dst = dstf;
        nbiter = 1;
    end
    ncoincs = dst.Struct.Setup.TotalCoinc;
    DetStruct = dst.Struct.Setup.Det;
    Detectors = [DetStruct.Name];
    CoincStruct = dst.Struct.Coinc;
    tag = CoincStruct.Det.Tag;
    amp = CoincStruct.Det.AmpMax;
    sig = CoincStruct.Det.Sigma;
    cordelays = 0;
    if cordelays == 1
        PlanStruct = CoincStruct.DelayCorrRecons.PlanRecons;
        SphStruct = CoincStruct.DelayCorrRecons.SphRecons;
    else
        PlanStruct = CoincStruct.PlanRecons.Radio;
        SphStruct = CoincStruct.SphRecons;
    end
    RunSetup = dst.Struct.Setup;
    X = [RunSetup.Det.X];
    Y = [RunSetup.Det.Y];
    Z = [RunSetup.Det.Z];
    mult = CoincStruct.Mult;
    timemat = CoincStruct.Det.UnixTime;
    times = max(timemat,[],2);  % seconds
    date_start = min(times);
    times = times-date_start;
    time = times./60; %minutes
    idp = CoincStruct.IdCoinc;
    thetap = PlanStruct.Theta;
    phip = PlanStruct.Phi;
    thetas = SphStruct.Theta;
    phis = SphStruct.Phi;  
    r = SphStruct.minDistSource;
    chi2s = SphStruct.Chi2Delay;
    slopes = SphStruct.SlopeDelay;
    chi2p = PlanStruct.Chi2Delay;
    slopep = PlanStruct.SlopeDelay;


    %% First selection
    disp(sprintf('%d coincs in total.',ncoincs))
    sel = find(mult>4);
    %disp(sprintf('With L>4: %d',length(sel)))
    sel = intersect(sel,find(chi2s<50 & abs(slopes-1)<1.1));
    %disp(sprintf('With valid spherical recons: %d',length(sel)))
    sel = intersect(sel, find(chi2p<50 & abs(slopep-1)<1.1));
    %disp(sprintf('With valid plan recons: %d',length(sel)))
    sel = intersect(sel, find(r>cutsettings.RCut));
    %disp(sprintf('With radius>%d m: %d',RCut,length(sel)))
    sel = intersect(sel, find(thetap<cutsettings.ThetaCut));
    %disp(sprintf('With ThetaPlan < %d deg: %d',ThetaCut,length(sel)))

    %% NS polar cut
    ewpolar = [148:155];
    [c,indEW] = intersect(Detectors,ewpolar);
    sel = intersect(sel,find(sum(tag(:,indEW),2)==0));
    %disp(sprintf('With no trigger on EW polar antennas: %d',length(sel)))

    %% Loop on candidates
    disp(sprintf('%d possible candidates to be checked.',length(sel)))
    for i=1:length(sel)
        if i/100==floor(i/100)
            disp(sprintf('%d/%d',i,length(sel)))
        end
        ind = sel(i);
        ind_tag_cand = find(tag(ind,:)>0);

        %% Barycenter
        bary_cand=[mean(X(ind_tag_cand)) mean(Y(ind_tag_cand)) mean(Z(ind_tag_cand))];
        barydist = zeros(1,length(ind_tag_cand));
        for m = 1:length(ind_tag_cand)
            barydist(m) = norm(bary_cand-[X(ind_tag_cand(m)) Y(ind_tag_cand(m)) Z(ind_tag_cand(m))]);
        end
        if mean(barydist)>cutsettings.BaryCut
          %disp(sprintf('Distance to barycenter > %d m. Skip candidate.',BaryCut))
          continue
        end
        
        %% Amplitude
        in = find(amp(ind,:)>0);
        calamp = amp(ind,in)./sig(ind,in);  % Calibrated amplitude (std dev normalisation)
        dets = Detectors(in);
        [maxi indmax] = max(calamp);
        mini = min(calamp);
        RatioAmp = maxi/mini;
        if RatioAmp<cutsettings.AmpRatioCut
            %disp(sprintf('Ratio MaxAmp/MinAmp<%3.1f. Skip candidate.',cutsettings.AmpRatioCut))
            continue
        end
        if abs(Y(in(indmax)))>100 | dets(indmax)>=148 | dets(indmax)==119 | dets(indmax)==120 
            %disp(sprintf('Max amplitude measured on external antenna (A%d). Skip candidate.',dets(indmax)))
            continue
        end
        
        %% Matrix of neighbourgh events       
        % Same ants
        TimeCut = [30 60 120];  %seconds
        com = zeros(3,5);
        frac = [0 0.3 0.5 0.66 0.8];
        for t = 1:length(TimeCut)
            timesel = find(abs(times-times(ind))<TimeCut(t));  %seconds
            common = zeros(1,length(timesel));
            for j = 1:length(timesel)
                toto=tag(timesel(j),:)+tag(ind,:);
                common(j) = length(find(toto==2))./mult(ind);
            end
            for k = length(frac):-1:1
                com(t,k) = length(find(common>=frac(k)))-1;
            end    
        end
        % Same dir
        azsel = find(chi2p<50 & abs(slopep-1)<1.1 & chi2s<50 & abs(slopes-1)<1.1 & r>500 & abs(phip-phip(ind))<cutsettings.PhiCut);
        DirTimeCut = [180 600 1200];  %seconds
        comdir = zeros(3,5);
        for t = 1:length(DirTimeCut)
            timesel2 = find(abs(times-times(ind))<DirTimeCut(t));  %seconds
            neighbourgs = intersect(timesel2,azsel);
            common = zeros(1,length(neighbourgs));
            for j = 1:length(neighbourgs)
                toto=tag(neighbourgs(j),:)+tag(ind,:);
                common(j) = length(find(toto==2))./mult(ind);
            end
            for k = length(frac):-1:1
                comdir(t,k) = length(find(common>=frac(k)))-1;
            end    
        end
        ncommon = com(find(TimeCut==cutsettings.TimeCut),find(frac==cutsettings.AntRatioCut));
        ncomdir = comdir(find(DirTimeCut/60==cutsettings.DirTimeCut),find(frac==cutsettings.AntRatioCut));
        if ncommon<=cutsettings.MaxAnt & ncomdir<=cutsettings.DirMaxAnt
            disp 'Candidate selected!'
            CandidateRun=[CandidateRun nrun];
            CandidateCoinc=[CandidateCoinc idp(ind)];
            CandidateTime=[CandidateTime max(timemat(ind,:))];
            CandidateThetaP=[CandidateThetaP thetap(ind)];
            CandidatePhiP=[CandidatePhiP phip(ind)];
            CandidateThetaS=[CandidateThetaS thetas(ind)];
            CandidatePhiS=[CandidatePhiS phis(ind)];
            CandidateRadius = [CandidateRadius r(ind)];
            CandidateNeighbourgs{end+1} = com;
            CandidateDirNeighbourgs{end+1} = comdir;
            CandidateAntennas{end+1} = Detectors(ind_tag_cand);
        end
    end
    meta = meta+1;
end
disp(sprintf('%d candidates selected in run %d.',length(find(CandidateRun==nrun)),nrun))

filename = 'Candidates_NS.mat';
save(filename,'cutsettings','CandidateRun','CandidateCoinc','CandidateTime','CandidateAntennas','CandidateNeighbourgs','CandidateDirNeighbourgs','CandidateThetaP','CandidatePhiP','CandidateThetaS','CandidatePhiS','CandidatePhiS','CandidateRadius')

fclose all;