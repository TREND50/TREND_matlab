function [] = SelectecCandidates(nrun,dstf)
% Select candidates in NS polar data
% OMH 11/06/2013

SharedGlobals;
DISPLAY = 0;
ThetaCut = 85;
RCut = 1000;  % % Mini radius [m]
BaryCut = 500; % Mean diustance to coinc barycenter [m]
TimeCut = 3; % Time window [minutes]
PhiCut = 10;  % [deg]
AntRatioCut = 2./3; % Max ratio of antennas in common with neighbouring events 
AmpRatioCut = 1.5;

%% Load list
listname = 'candidatesNS.mat';
if fopen(listname)>0
    list = open(listname);
    candrunid = list.runid;
    candcoincid = list.coincid;
else
    candrunid = [];
    candcoincid = [];
end
candcoincid(find(candrunid==nrun)) = [];
candrunid(find(candrunid==nrun)) = [];

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
        dstname = [DST_PATH sprintf(dst_filename,nrun,meta)];
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
    time = max(CoincStruct.Det.UnixTime,[],2)/60;  % minutes
    date_start = min(time);
    time = time-date_start;
    idp = CoincStruct.IdCoinc;
    thetap = PlanStruct.Theta;
    phip = PlanStruct.Phi;
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
    sel = intersect(sel, find(r>RCut));
    %disp(sprintf('With radius>%d m: %d',RCut,length(sel)))
    sel = intersect(sel, find(thetap<ThetaCut));
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
            disp(sprintf('Checking coinc %d: t = %3.1f mins, thetap = %3.1f deg, phip = %3.1f deg...',idp(ind),time(ind),thetap(ind),phip(ind)))  
            disp 'Detectors:'
            Detectors(ind_tag_cand)
            disp(sprintf('Distance to barycenter = %3.1f m. ',mean(barydist)))
        end
        if mean(barydist)>BaryCut
          %disp(sprintf('Distance to barycenter > %d m. Skip candidate.',BaryCut))
          continue
        end
        
        % Amplitude
        calamp = amp(ind,:)./sig(ind,:);  % Calibrated amplitude (std dev normalisation)
        [maxi indmax] = max(calamp);
        mini = min(calamp);
        RatioAmp = maxi/mini;
        if RatioAmp<AmpRatioCut
            disp(sprintf('Ratio MaxAmp/MinAmp<%3.1f. Skip candidate.',AmpRatioCut))
            continue
        end
        if abs(Y(indmax))>100 | Detectors(indmax)>=156 | Detectors(indmax)==119 | Detectors(indmax)==120 
            disp(sprintf('Max amplitude measured on external antenna (A%d). Skip candidate.',Detectors(indmax)))
            continue
        end

        % Neighbourgs
        timesel = find(abs(time-time(ind))<TimeCut);
        azsel = find(r>300 & abs(phip-phip(ind))<PhiCut);
        neighbourgs = intersect(timesel,azsel);
        common = zeros(1,length(neighbourgs));
        for j = 1:length(neighbourgs)
            toto=tag(neighbourgs(j),:)+tag(ind,:);
            common(j) = length(find(toto==2))./mult(ind);
        end
        ncommon = length(find(common>=AntRatioCut));
        if ncommon<2
            if DISPLAY
                disp(sprintf('Coincs in +- %d mins around candidate: %d.',TimeCut,length(timesel)-1))
                disp(sprintf('And in +- %d deg in azimuth: %d.',PhiCut,length(neighbourgs)-1))
                disp(sprintf('And with at least %3.2f of antennas in common with candidate: %d.',AntRatioCut,ncommon-1))
                disp(sprintf('Ratio MaxAmp/MinAmp=%3.1f. OK.',RatioAmp))
                disp(sprintf('Max amplitude measured on  antenna %d. OK.',Detectors(indmax)))
                %
                AnaCoinc(nrun,idp(ind),dst.Struct)
                pause
                close all
            end
            candrunid(end+1) = nrun;
            candcoincid(end+1) = idp(ind);
        else
            if DISPLAY
              disp(sprintf('%d neighbour coincs with similar pattern. Skip candidate.',ncommon-1))
            end
        end
    end
    meta = meta+1;
end
disp(sprintf('%d candidates selected in run %d.',length(find(candrunid==nrun)),nrun))
runid = candrunid;
coincid = candcoincid;
save('candidatesNS.mat','runid','coincid');
