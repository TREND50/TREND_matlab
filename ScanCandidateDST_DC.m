function [thpel phpel mult radel rampel run coinc val] = ScanCandidateDST_DC(E,PeriodId,down)
% Adjust cuts on neighbourgs in Candidate DST
% OMH 21/06/2013
% down = 1: copy selection to CandidateTag.mat  [DANGEROUS!]
% down = 0: only displays all selected events
% down = -1: dispay selected==1 AND tag==1 events (Double selection)
% Adapted for DC
% SL 01/04/2017

SharedGlobals;

if ~exist('down')  
  down = 0;
end

%E='3e18'

DCfile=load([CAND_PATH '/' E '/log_candanalysis.txt'], 'r' );
runs=DCfile(:,4);
simucoincs=DCfile(:,13);
tag=DCfile(:,25);

dstname = [CAND_PATH E sprintf('/Candidates_Period%d_102014.mat',PeriodId)];
if fopen(dstname)<0
    disp(sprintf('File %s not found. Abort.',dstname))
    thpel = [];
    phpel = [];
    mult = [];
    RunId = [];
    CoincId = [];
    return
end
c = load(dstname);
run = c.CandidateRun;
ncand = length(run);
coinc = c.CandidateCoinc;
ants = c.CandidateAntennas;
nei = c.CandidateNeighbourgs;
dirnei = c.CandidateDirNeighbourgs;
thetap = c.CandidateThetaP;
phip = c.CandidatePhiP;
thetas = c.CandidateThetaS;
phis = c.CandidatePhiS;
chi2p = c.CandidateChi2P;
bad = c.CandidateBadAnt;
radius = c.CandidateRadius;
ramp = c.CandidateRatioAmp;
valid = zeros(1,ncand);
maxmult = 0;
n1 = 0;
n2 = 0;
n3 = 0;
n4 = 0;
n5 = 0;
n6 = 0;
n7 = 0;
n8 = 0;

DCfilename = [CAND_PATH E '/log_candanalysis.txt']; 
DCfile=load(DCfilename, 'r' );
fid = fopen([CAND_PATH E '/log_candanalysisend.txt'], 'a' ); %SL
for i = 1:ncand
    %disp(sprintf('R%dC%d',run(i),coinc(i)))
    if bad(i)>0
        %disp 'Skip 1.'
        n1 = n1+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'0 0 0 0 0 0 0 0\n')
        continue
    end
    
    if chi2p(i)>30
        %disp 'Skip 2.'
        n2 = n2+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 0 0 0 0 0 0 0\n')
        continue
    end
        
    if radius(i)<3000
        %disp 'Skip 6.'
        n6 = n6+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 0 0 0 0 0 0\n')
        continue
    end
    
    if ramp(i)<1
        n7 = n7+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 1 0 0 0 0 0\n')
        continue
    end
    
    %if thetap(i)>80 | thetap(i)<0
    if thetap(i)>81 | thetap(i)<0      %SL test
        n3 = n3+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 1 1 0 0 0 0\n')
        continue
    end
    
    mdirnei = dirnei{i};
%    if mdirnei(4,1)>0  % Hardest
%    if mdirnei(1,5)>0  % Softest
    if mdirnei(3,2)>0  %Best 10mn + 33%
%    if mdirnei(3,1)>1  %EW
%    if mdirnei(3,2)>1  %NS
        %disp 'Skip 4.'
        n4 = n4+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 1 1 1 0 0 0\n')
        continue
    end
    
    mnei = nei{i};
%    if mnei(4,1)>0  % Hardest 
%    if mnei(1,5)>0  %Softest
    if mnei(1,4)>0  %Best: 30s+66%
%    if mnei(1,4)>0  %EW   30s+66%
%    if mnei(1,3)>0  %NS
        %disp 'Skip 5.'
        n5 = n5+1;
        thiscoinc=coinc(i);
        thisind=find(simucoincs==thiscoinc);
        thisline=DCfile(thisind,:);
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 1 1 1 1 0 0\n')
        continue
    end
    
    if PeriodId>=9
        mants = ants{i};
        badpat = sum(ismember(mants,[129 131 135 136 137 138]));
        if badpat>=2
            disp 'Skip 8.'
            n8 = n8+1;
            thiscoinc=coinc(i);
            thisind=find(simucoincs==thiscoinc);
            thisline=DCfile(thisind,:);
            fprintf(fid,'%d ',thisline);
            fprintf(fid,'1 1 1 1 1 1 1 0\n')
            continue
        end
        ewpol = sum(ismember(mants,[148:153]));  % Antennas measuring EW polar
        if ewpol>=1
            disp 'Skip 9.'
            continue
        end
    end    
    
    thiscoinc=coinc(i);
    thisind=find(simucoincs==thiscoinc);
    thisline=DCfile(thisind,:);
    if size(thisline,1) == 0
        disp(sprintf('Error! Candidate Coinc Id %d not found in %s',thiscoinc,DCfilename))
        pause
    else
        fprintf(fid,'%d ',thisline);
        fprintf(fid,'1 1 1 1 1 1 1 1\n')
        disp 'Selected.'
        valid(i) = 1;
    end
end
sel = find(valid==1);
val = [ncand n1 n2 n6 n7 n3 n4 n5 n8 length(sel)]

disp(sprintf('%d candidates in DST, %d selected.',ncand,length(sel)))

figure(1)
subplot(2,1,1)
hist(thetap(sel),20)
xlim([0 90])
subplot(2,1,2)
hist(phip(sel),40)
xlim([0 360])

figure(2)
PrepareSkyPlot(2);
polar( phip(sel )*DEG2RAD(1), thetap(sel), 'go' );  

thpel = thetap(sel);
phpel = phip(sel);
radel = radius(sel);
rampel = ramp(sel);
run = run(sel);
coinc = coinc(sel);


[run' coinc']
