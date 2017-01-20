function [] = CandidateAnalysis_NS(nrun,dstf)
% Select candidates in NS polar data
% Writes it to dst.
% OMH 11/06/2013

SharedGlobals;
DISPLAY = 0;
cutsettings.ThetaCut = 85;
cutsettings.RCut = 1000;  % % Mini radius [m]
cutsettings.BaryCut = 500; % Mean diustance to coinc barycenter [m]
cutsettings.TimeCut = 30; % Time window [seconds]
cutsettings.DirTimeCut = 10; % Time window for same direction [minutes]
cutsettings.MaxAnt = 10;
cutsettings.DirMaxAnt = 2;
cutsettings.PhiCut = 10;  % Same direction [deg]
cutsettings.AntRatioCut = 0.66; % Max ratio of antennas in common with neighbouring events 
cutsettings.AmpRatioCut = 1.5;


%% Load list
filename = 'Candidates_NS.mat';
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
    %
    %Clear data for this run
    thisRun = find(CandidateRun==nrun);
    CandidateThetaS(thisRun) = [];
    CandidatePhiS(thisRun) = [];
    CandidateRadius(thisRun) = [];
    CandidateAntennas(thisRun) = [];
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
end


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

        % Barycenter
        bary_cand=[mean(X(ind_tag_cand)) mean(Y(ind_tag_cand)) mean(Z(ind_tag_cand))];
        barydist = zeros(1,length(ind_tag_cand));
        for m = 1:length(ind_tag_cand)
            barydist(m) = norm(bary_cand-[X(ind_tag_cand(m)) Y(ind_tag_cand(m)) Z(ind_tag_cand(m))]);
        end
        if DISPLAY 
            disp(sprintf('\n**Checking coinc %d: \nt = %3.1f mins, thetap = %3.1f deg, phip = %3.1f deg...',idp(ind),time(ind),thetap(ind),phip(ind)))  
            disp 'Detectors:'
            Detectors(ind_tag_cand)
            disp(sprintf('Distance to barycenter = %3.1f m. ',mean(barydist)))
        end
        if mean(barydist)>cutsettings.BaryCut
          %disp(sprintf('Distance to barycenter > %d m. Skip candidate.',BaryCut))
          continue
        end
        
        % Amplitude
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
        
        % Neighbourgs
        % Same ants
        timesel = find(abs(times-times(ind))<cutsettings.TimeCut);
        common = zeros(1,length(timesel));
        for j = 1:length(timesel)
            toto=tag(timesel(j),:)+tag(ind,:);
            common(j) = length(find(toto==2))./mult(ind);
        end
        ncommon = length(find(common>=cutsettings.AntRatioCut))-1;
        % Same dir
        azsel = find(r>500 & abs(phip-phip(ind))<cutsettings.PhiCut);
        timesel2 = find(abs(time-time(ind))<cutsettings.DirTimeCut);
        neighbourgs = intersect(timesel2,azsel);
        comdir = zeros(1,length(neighbourgs));
        for j = 1:length(neighbourgs)
            toto=tag(neighbourgs(j),:)+tag(ind,:);
            comdir(j) = length(find(toto==2))./mult(ind);
        end
        ncomdir = length(find(comdir>=cutsettings.AntRatioCut))-1;
        if ncommon<=cutsettings.MaxAnt & ncomdir<=cutsettings.DirMaxAnt
            if DISPLAY
                disp(sprintf('Coincs in +- %d seconds around candidate: %d.',cutsettings.TimeCut,length(timesel)))
                disp(sprintf('And with at least %3.2f pc of antennas in common with candidate: %d.',cutsettings.AntRatioCut*100,ncommon))
                disp(sprintf('Coincs in +- %d minutes around candidate in +- %d deg in azimuth: %d.',cutsettings.PhiCut,cutsettings.DirTimeCut,length(neighbourgs)))
                disp(sprintf('And with at least %3.2f pc of antennas in common with candidate: %d.',cutsettings.AntRatioCut*100,ncomdir))
                disp(sprintf('Ratio MaxAmp/MinAmp=%3.1f. OK.',RatioAmp))
                disp(sprintf('Max amplitude measured on  antenna %d. OK.',Detectors(indmax)))
                %
                %AnaCoinc(nrun,idp(ind),dst.Struct)
                pause
                close all
            end
            disp 'Candidate selected!'
            %pause
            CandidateRun=[CandidateRun nrun];
            CandidateCoinc=[CandidateCoinc idp(ind)];
            CandidateTime=[CandidateTime max(timemat(ind,:))];
            CandidateThetaP=[CandidateThetaP thetap(ind)];
            CandidatePhiP=[CandidatePhiP phip(ind)];
            CandidateThetaS=[CandidateThetaS thetas(ind)];
            CandidatePhiS=[CandidatePhiS phis(ind)];
            CandidateRadius = [CandidateRadius r(ind)];
            [CandidateAntennas{end+1}] = Detectors(ind_tag_cand);
        else
            if DISPLAY
              disp(sprintf('%d neighbour coincs with similar pattern, %d with same dir. Skip candidate.',ncommon,ncomdir))
            end
        end
    end
    meta = meta+1;
end
disp(sprintf('%d candidates selected in run %d.',length(find(CandidateRun==nrun)),nrun))

filename = 'Candidates_NS.mat';
save(filename,'cutsettings','CandidateRun','CandidateCoinc','CandidateTime','CandidateAntennas','CandidateThetaP','CandidatePhiP','CandidateThetaS','CandidatePhiS','CandidatePhiS','CandidateRadius')

fclose all;