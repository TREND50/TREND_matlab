function [] = CandidateAnalysis_DC(E,nrun,dstf)
% Select candidates v2.
% Writes it to dst.
% OMH 18/06/2013
% Adapted to DataChallenge analysis
% SL 01/04/2017

SharedGlobals;
if DC==0
    disp (sprintf('Error! DC = %d! Update SharedGlobals.m Aborting.',DC))
    return
end
%E='3e18' %SL
DISPLAY=0 %SL
if nrun == 3202  
  disp 'Skip R3202, too big'
  return
end
isSimu = 0;
cutsettings.ThetaCut = 81; %SL test
cutsettings.RCut = 500;  % % Mini radius [m]
cutsettings.Chi2sCut = 30;
cutsettings.Chi2pCut = 30;
%cutsettings.Chi2pCut = 50;
cutsettings.BaryCut = 500; % Mean distance to coinc barycenter [m]
cutsettings.AmpRatioCut = 1.0; %ORIGINAL = 1.5
cutsettings.TimeCut = 30; % Time window [seconds]
cutsettings.DirTimeCut = 3; % Time window for same direction [minutes]
cutsettings.MaxAnt = 10;
cutsettings.PhiCut = 10;  % Same direction [deg]
cutsettings.DirMaxAnt = 3;
cutsettings.AntRatioCut = 0.66; % Max ratio of antennas in common with neighbouring events 
cutsettings.TimeVector = [30 60 120 150];  %seconds
cutsettings.DirTimeVector = [180 600 1200 1800];  %seconds
cutsettings.Frac = [0 0.3 0.5 0.66 0.8]; % Fraction of antennas in common
frac = cutsettings.Frac;
alivecut = 100;  

%% Read simu recover %SL
DCfile=load([DST_PATH E '/log_recoverall.txt'], 'r' );
runs=DCfile(:,4);
simucoincs=DCfile(:,13);
trigscore=DCfile(:,8);
simucoinc=simucoincs(simucoincs~=0 & runs==nrun & trigscore>=4);

%% Load valid Candidates DST
periods = zeros(10,2);
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

periodID = find(periods(:,1)<=nrun & periods(:,2)>=nrun);

if size(periodID,1)==0
    disp(sprintf('Run %d does not belong to an identified period! Abort.',nrun))
    return
end

candname = [CAND_PATH E '/' sprintf('Candidates_Period%d_102014.mat',periodID)];
disp(sprintf('Now loading DST %s...',candname))
if fopen(candname)>0
    % load dst
    c = open(candname);
    CandidateRun = c.CandidateRun;
    CandidateCoinc = c.CandidateCoinc;
    CandidateTime = c.CandidateTime;
    CandidateThetaP = c.CandidateThetaP;
    CandidatePhiP = c.CandidatePhiP;
    CandidateChi2P = c.CandidateChi2P;    
    CandidateThetaS = c.CandidateThetaS;
    CandidatePhiS = c.CandidatePhiS;
    CandidateRadius = c.CandidateRadius;
    CandidateAntennas = c.CandidateAntennas;
    CandidateMaxAnt = c.CandidateMaxAnt;
    CandidateBadAnt = c.CandidateBadAnt;   
    CandidateNeighbourgs = c.CandidateNeighbourgs;
    CandidateDirNeighbourgs = c.CandidateDirNeighbourgs;
    CandidateRatioAmp = c.CandidateRatioAmp;
    %Clear data for this run
    thisRun = find(CandidateRun==nrun);
   if size(thisRun,2)
       disp(sprintf('Run %d already present, skip.',nrun))
       fclose all;
       return
   end
else
    disp(sprintf('No DST %s...',candname))
    CandidateRun=[];
    CandidateCoinc=[];
    CandidateTime=[];
    CandidateThetaP=[];
    CandidatePhiP=[];
    CandidateChi2P=[];
    CandidateThetaS=[];
    CandidatePhiS=[];
    CandidateRadius=[];
    CandidateMaxAnt=[];
    CandidateBadAnt=[];
    CandidateAntennas={};
    CandidateNeighbourgs={};
    CandidateDirNeighbourgs={};
    CandidateRatioAmp = [];   
end
disp 'Done.'

%% Get number of sub dsts
stopflag=0;
nbiter=1;
while stopflag==0
    filename = [DST_PATH E sprintf('/dst%d_%d.mat',nrun,nbiter)];
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

DCfile=load([DST_PATH E '/log_recoverall.txt'], 'r' );
fid = fopen([CAND_PATH E '/log_candanalysis.txt'], 'a+' ); %SL

%% Loop on sub dsts
meta = 1;
while meta<=nbiter

    %% Load dst
    if ~exist('dstf')
        dstname = [DST_PATH E sprintf('/dst%d_%d.mat',nrun,meta)];
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
    evt = [DetStruct.Evt];
    CoincStruct = dst.Struct.Coinc;
    tag = CoincStruct.Det.Tag;
    amp = CoincStruct.Det.AmpMax;
    sig = CoincStruct.Det.Sigma;
    if isSimu
        stat = ones(ncoincs,length(Detectors));
    else
        stat = CoincStruct.Det.Status;
    end
    PlanStruct = CoincStruct.PlanRecons.Radio;
    SphStruct = CoincStruct.SphRecons;
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
    xs = SphStruct.X0;
    ys = SphStruct.Y0;
    zs = SphStruct.Z0;
    phip = PlanStruct.Phi;
    r = SphStruct.minDistSource;
    chi2s = SphStruct.Chi2Delay;
    slopes = SphStruct.SlopeDelay;
    chi2p = PlanStruct.Chi2Delay;
    slopep = PlanStruct.SlopeDelay;
    aDetectors = Detectors(evt>alivecut); % "Alive" detectors
    aX = X(evt>alivecut);
    aY = Y(evt>alivecut);
    aZ = Z(evt>alivecut);
    aevt = evt(evt>alivecut);
    cen = aDetectors(find(abs(aY)<100));
    cen(find(cen==119)) = [];
    cen(find(cen==120)) = [];
    %cen(find(cen==137)) = [];
    cen(end+1)=138;
    cen(end+1)=140;
    [a icen] = intersect(aDetectors,cen);
    neid = {};
    for i = 1:length(icen)
        neid(i,1) = {aDetectors(icen(i))};
        delta = [aX'-aX(icen(i)) aY'-aY(icen(i)) aZ'-aZ(icen(i))]; 
        dist = sqrt(sum(delta.*delta,2));
        [dist idist] = sort(dist);
        nmax = min(22,length(aDetectors));
        nei = aDetectors(idist(2:nmax));
        philoc = mod(atan2(-(aX(idist(2:nmax))-aX(icen(i))),(aY(idist(2:nmax))-aY(icen(i))))*180/pi,360);        
        for j = 1:4
            if j==1 | j==3  % Select detectors in North & South quadrants
                sel = find(abs(philoc-90*(j-1))<45);  
            else  % East & West
                 sel = find(abs(philoc-90*(j-1))<45);
            end
            neid(i,j+1) = {nei(sel)};
        end
    end
    
    %% First selection
    disp(sprintf('%d coincs in total.',ncoincs))
    disp(sprintf('%d coincs are simu coinc.',length(simucoinc)))
    [a indsim b ] = intersect(idp,simucoinc);
    [a inddata b] = setxor(idp,simucoinc);
    sel = find(mult>4);
    simusel= intersect(sel, indsim); %SL
    sel=intersect(sel,inddata);
    disp(sprintf('Data With L>4: %d',length(sel)))
    disp(sprintf('Simu With L>4: %d',length(simusel)))%SL
    res = [length(inddata) length(sel)];
    simures=[length(indsim) length(simusel)]; %SL
    
    for s=1:length(indsim) %SL
        inter=intersect(indsim(s),simusel);
        if isempty(inter)
            thiscoinc = idp(indsim(s));
            thisind = find(simucoincs==thiscoinc);
            thisline = DCfile(thisind,:);
            fprintf(fid,'%d ',thisline)
            fprintf(fid,'0 0 0 0 0 0 0 0 0 \n')
        end
    end
    
    for s=1:length(simusel) %SL
        inter=intersect(simusel(s),find(r>cutsettings.RCut));
        if isempty(inter)
            thiscoinc = idp(simusel(s));
            thisind = find(simucoincs==thiscoinc);
            thisline = DCfile(thisind,:);
            fprintf(fid,'%d ',thisline)
            fprintf(fid,'1 0 0 0 0 0 0 0 0\n')  % Setting field 1 to 1 <=> validating previous cut (Mult). Field 2 (Rcut) 1st field at 0 
        end
    end
    sel = intersect(sel, find(r>cutsettings.RCut));
    simusel=intersect(simusel, find(r>cutsettings.RCut));%SL
    disp(sprintf('Data With radius>%d m: %d',cutsettings.RCut,length(sel)))
    disp(sprintf('Simu With radius>%d m: %d',cutsettings.RCut,length(simusel)))%SL
    res = [res length(sel)];
    simures=[simures length(simusel)]; %SL    
    for s=1:length(simusel) %SL
        inter=intersect(simusel(s),find(chi2p<cutsettings.Chi2pCut));
        if isempty(inter)
            thiscoinc=idp(simusel(s));
            thisind=find(simucoincs==thiscoinc);
            thisline=DCfile(thisind,:);
            fprintf(fid,'%d ',thisline)
            fprintf(fid,'1 1 0 0 0 0 0 0 0\n')
        end
    end
    sel = intersect(sel, find(chi2p<cutsettings.Chi2pCut));
    simusel = intersect(simusel, find(chi2p<cutsettings.Chi2pCut)); %SL
    res = [res length(sel)];
    simures = [simures length(simusel)];%SL
    disp(sprintf('Data With valid plan recons: %d',length(sel)))
    disp(sprintf('Simu With valid plan recons: %d',length(simusel)))%SL  
    

    for s=1:length(simusel) %SL
        inter=intersect(simusel(s),find(thetap<cutsettings.ThetaCut));
        if isempty(inter)
            thiscoinc=idp(simusel(s));
            thisind=find(simucoincs==thiscoinc);
            thisline=DCfile(thisind,:);
            fprintf(fid,'%d ',thisline)
            fprintf(fid,'1 1 1 0 0 0 0 0 0 \n')
        end
    end
    sel = intersect(sel, find(thetap<cutsettings.ThetaCut));
    simusel = intersect(simusel, find(thetap<cutsettings.ThetaCut));%SL
    disp(sprintf('Data With ThetaPlan < %d deg: %d',cutsettings.ThetaCut,length(sel)));
    disp(sprintf('Simu With ThetaPlan < %d deg: %d',cutsettings.ThetaCut,length(simusel)));%SL
    res = [res length(sel)];
    simures = [simures length(simusel)];%SL
         
    %% Loop on candidates
    if DC %SL
        res=simures
        sel=simusel
    end
    disp(sprintf('%d possible candidates to be checked.',length(sel)))
    %pause
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n4 = 0;
    n5 = 0;
    for i=1:length(sel)
        if i/100==floor(i/100)
            disp(sprintf('%d/%d',i,length(sel)))
        end
        ind = sel(i);
        in = find(tag(ind,:)>0);
        out = find(tag(ind,:)==0);
        
        %% Barycenter
        bary_cand=[mean(X(in)) mean(Y(in)) mean(Z(in))];
        [t,p,rc] = cart2sph(xs(ind)-bary_cand(1),ys(ind)-bary_cand(2),zs(ind)-bary_cand(3));
        phisc = t*180/pi-90;
        phisc = mod(phisc,360);
        thsc = 90-p*180/pi;
        barydist = zeros(1,length(in));
        for m = 1:length(in)
            barydist(m) = norm(bary_cand-[X(in(m)) Y(in(m)) Z(in(m))]);
        end
        if mean(barydist)>cutsettings.BaryCut
          %disp(sprintf('Coinc %d: distance to barycenter > %d m. Skip candidate.',idp(ind),cutsettings.BaryCut))
          n1 = n1+1;
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline)
              fprintf(fid,'1 1 1 1 0 0 0 0 0\n')
            end
          continue
        end
        
        %% Amplitude
        calamp = amp(ind,in)./sig(ind,in);  % Calibrated amplitude (std dev normalisation)
        dets = Detectors(in);
        [maxi indmax] = max(calamp);
        mini = min(calamp);
        RatioAmp = maxi/mini;
        
        %% Status
        status = stat(ind,:);
        thisStatus = zeros(length(Detectors),10);
        thisStatus(:,1) = Detectors';    
        for j = 1:length(Detectors)
          st = dec2bin(status(j));
          stv = int16(sscanf(st,'%1d'))'; %LSB last
          stv=fliplr(stv); %LSB first
          thisStatus(j,2:length(stv)+1) = stv;
        end
        bads = (thisStatus(:,3)==1);
        if sum(bads)>1
            %disp(sprintf('Coinc %d: bad signals on %d antennas. Skip candidate.',idp(ind),sum(bads)))
            n2 = n2+1;
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline);
              fprintf(fid,'1 1 1 1 1 0 0 0 0\n')
            end
            continue
        end
        
        %% Pattern at ground
        [a aout] = intersect(aDetectors,Detectors(out));
        icenout = intersect(aout,icen); % indexes of non-triggered inner antennas
        [a ineid] = intersect([neid{:,1}],aDetectors(icenout));
        skip = 0;
        j = 1;
        while j<=length(icenout) & skip<2 % Loop on non-triggered inner antennas
          if aevt(icenout(j))<alivecut % Dead antenna: skip it 
              j = j+1;
              disp 'problem! Should not be here'
              pause
          end
          if size(neid{ineid(j),3},2) == 0  % No detector in the left dir
              lin = 1;   
          else
              lin = length(intersect(Detectors(in),neid{ineid(j),3}));  % Triggered detectors in left direction
          end
          if size(neid{ineid(j),5},2) == 0
              rin = 1;
          else
              rin = length(intersect(Detectors(in),neid{ineid(j),5}));  % Same in Right
          end
          if size(neid{ineid(j),2},2) == 0
              uin = 1;
          else
              uin = length(intersect(Detectors(in),neid{ineid(j),2}));  % Up
          end          
          if size(neid{ineid(j),4},2) == 0
              din = 1;
          else
              din = length(intersect(Detectors(in),neid{ineid(j),4}));  % Down
          end
          %[neid{ineid(j),1} lin rin uin din]
          %neid{ineid(j),2}
          if (lin>0 & rin>0) & (uin>0 & din>0)          
              %disp 'Hole!'
              %disp(sprintf('Wrong ground pattern: antenna %d did not trigger. ',Detectors(icenout(j))))
              skip = skip+1;
              %break
          end
          j = j+1;
        end
        %
        if DISPLAY == 1
            figure(1)
            hold off
            plot(X,Y,'x')
            hold on
            plot(aX,aY,'^','MarkerFaceColor','y')
            grid on
            plot(X(icenout),Y(icenout),'^','MarkerFaceColor','k')
            plot(X(in),Y(in),'^','MarkerFaceColor','g','MarkerSize',12)
            xmin = -100;
            xmax = 3000;
            ymin = -500;
            ymax = 800;
            axis([xmin,xmax,ymin,ymax])
            axis equal
            for k = 1:length( Detectors )
                text( X(k)+20, Y(k), num2str( Detectors(k) ), 'FontSize', 12, 'FontWeight', 'bold','Color','k' );
            end
        end
        if skip > 1           
            n3 = n3+1;
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline)
              fprintf(fid,'1 1 1 1 1 1 0 0 0 \n')
             end         
             continue
        end
        

        %% Matrix of neighbourgh events
        % Same ants
        com = zeros(3,5);
        TimeCut = cutsettings.TimeVector;
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
        if size(r,1)>1  % Radius defined as a  line vector instead of column for some runs
            r = r';
        end
        azsel = find(chi2p<cutsettings.Chi2pCut & abs(slopep-1)<1.1 & chi2s<cutsettings.Chi2sCut & abs(slopes-1)<1.1 & r>500 & abs(phip-phip(ind))<cutsettings.PhiCut);
        
        comdir = zeros(3,5);
        DirTimeCut = cutsettings.DirTimeVector;
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
        
        if (ncommon<=cutsettings.MaxAnt & ncomdir<=cutsettings.DirMaxAnt)
            disp(sprintf('Coinc %d: candidate selected!',idp(ind)))
            CandidateRun=[CandidateRun nrun];
            CandidateCoinc=[CandidateCoinc idp(ind)];
            CandidateTime=[CandidateTime max(timemat(ind,:))];
            CandidateThetaP=[CandidateThetaP thetap(ind)];
            CandidatePhiP=[CandidatePhiP phip(ind)];
            CandidateChi2P=[CandidateChi2P chi2p(ind)];
            CandidateThetaS=[CandidateThetaS thsc];
            CandidatePhiS=[CandidatePhiS phisc];
            CandidateMaxAnt=[CandidateMaxAnt dets(indmax)];
            CandidateBadAnt=[CandidateBadAnt sum(bads)];
            CandidateRadius = [CandidateRadius r(ind)];
            CandidateNeighbourgs{end+1} = com;
            CandidateRatioAmp=[CandidateRatioAmp RatioAmp];
            CandidateDirNeighbourgs{end+1} = comdir;
            CandidateAntennas{end+1} = Detectors(in);
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline);
              fprintf(fid,'1 1 1 1 1 1 1 1 1 \n')
            end   
         elseif  ncomdir>cutsettings.DirMaxAnt   
            n4 = n4+1;
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline)
              fprintf(fid,'1 1 1 1 1 1 1 0 0\n')
            end   
         elseif ncommon>cutsettings.MaxAnt   
            n5 = n5+1;
            if DC
              thiscoinc=idp(ind);
              thisind=find(simucoincs==thiscoinc);
              thisline=DCfile(thisind,:);
              fprintf(fid,'%d ',thisline)
              fprintf(fid,'1 1 1 1 1 1 1 1 0 \n')
            end   
        end
        %fclose(fid)
    end
    res = [nrun meta res n1 n2 n3 n4 n5 length(find(CandidateRun==nrun))]
    fid2 = fopen([DST_PATH '/' E '/' sprintf('selection_%d.txt',periodID)], 'a+' );
    fprintf( fid2, '%6d ', res(1));   % Run
    fprintf( fid2, '%6d ', res(2));   % MetaRun
    fprintf( fid2, '%6d ', res(3));   % Nini
    fprintf( fid2, '%6d ', res(4));   % Mult
    fprintf( fid2, '%6d ', res(5));   % Radius
    fprintf( fid2, '%6d ', res(6));   % ValidPlan
    fprintf( fid2, '%6d ', res(7));   % Theta
    fprintf( fid2, '%6d ', res(8));   % n1
    fprintf( fid2, '%6d ', res(9));   % n2
    fprintf( fid2, '%6d ', res(10)); % n3
    fprintf( fid2, '%6d ', res(11)); % n4
    fprintf( fid2, '%6d ', res(12));   % n5
    fprintf( fid2, '%6d ', res(13));   % nCands
    fprintf( fid2, '\n' );
    fclose(fid2);
    meta = meta+1;
end
fclose(fid);


disp(sprintf('%d candidates selected in run %d.',length(find(CandidateRun==nrun)),nrun))

if length(find(CandidateRun==nrun))>0  % Candidates have been found
    [CandidateRun ind] = sortrows(CandidateRun');
    if size(CandidateRun,1)>1
        CandidateRun = CandidateRun';
    end
    CandidateCoinc = CandidateCoinc(ind);
    CandidateTime = CandidateTime(ind);
    CandidateThetaP = CandidateThetaP(ind);
    CandidatePhiP = CandidatePhiP(ind);
    CandidateChi2P = CandidateChi2P(ind);
    CandidateThetaS = CandidateThetaS(ind);
    CandidatePhiS = CandidatePhiS(ind);
    CandidateMaxAnt = CandidateMaxAnt(ind);
    CandidateBadAnt = CandidateBadAnt(ind);
    CandidateRadius = CandidateRadius(ind);
    CandidateNeighbourgs = CandidateNeighbourgs(ind);
    CandidateRatioAmp = CandidateRatioAmp(ind);
    CandidateDirNeighbourgs = CandidateDirNeighbourgs(ind);
    CandidateAntennas = CandidateAntennas(ind);
    save(candname,'cutsettings','CandidateRun','CandidateCoinc','CandidateTime','CandidateRatioAmp','CandidateMaxAnt','CandidateBadAnt','CandidateAntennas','CandidateNeighbourgs','CandidateDirNeighbourgs','CandidateThetaP','CandidatePhiP','CandidateChi2P','CandidateThetaS','CandidatePhiS','CandidatePhiS','CandidateRadius')
end

fclose all;

