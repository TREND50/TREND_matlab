function [thpel phpel mult radel rampel run coinc val] = ScanCandidateDST(PeriodId,down)
% Adjust cuts on neighbourgs in Candidate DST
% OMH 21/06/2013
% down = 1: copy selection to CandidateTag.mat  [DANGEROUS!]
% down = 0: only displays all selected events
% down = -1: dispay selected==1 AND tag==1 events (Double selection)

if ~exist('down')  
  down = 0;
end

SharedGlobals;
%close all;

dstname = sprintf('Candidates_Period%d.mat',PeriodId);
dstname = [CAND_PATH dstname]; 
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

for i = 1:ncand
    %disp(sprintf('R%dC%d',run(i),coinc(i)))
    if bad(i)>0
        %disp 'Skip 1.'
        n1 = n1+1;
        continue
    end
    
    if chi2p(i)>30
        %disp 'Skip 2.'
        n2 = n2+1;
        continue
    end
        
    if radius(i)<3000
        %disp 'Skip 6.'
        n6 = n6+1;
        continue
    end
    
    if ramp(i)<1
        %disp 'Skip 7.'
        n7 = n7+1;
        continue
    end
    
    if thetap(i)>80 | thetap(i)<0
        %| abs(thetas(i)-thetap(i))>10 | abs(phis(i)-phip(i))>3 
        %disp 'Skip 3'
        n3 = n3+1;
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
        continue
    end
    
    if PeriodId>=9
        mants = ants{i};
        badpat = sum(ismember(mants,[129 131 135 136 137 138]));
        if badpat>=2
            disp 'Skip 8.'
            n8 = n8+1;
            continue
        end
       ewpol = sum(ismember(mants,[148:153]));  % Antennas measuring EW polar
        if ewpol>=1
            disp 'Skip 9.'
            continue
        end
    end    
    
    disp 'Selected.'
%     if length(ants{i}) > maxmult
%         maxmult = length(ants{i})
%         ants{i}
%         pause
%     end
    valid(i) = 1;
end
sel = find(valid==1);
val = [ncand n1 n2 n6 n7 n3 n4 n5 n8 length(sel)]
% ncand
% 
%pause

%% Load CandidateTag.mat
if fopen('CandidateTag.mat')>0
	t = load('CandidateTag.mat');
	runt = t.RunId;
	coinct = t.CoincId;
	tag = t.Tag;
else
	runt = [];
	coinct = [];
	tag = [];
end
sel2 = [];

if down == 1
    down
    disp 'Now validating all selected candidates (ie modification of CandidateTag.mat), OK?'
    pause
end

for i = 1:length(sel)
    ind = find(runt==run(sel(i)) & coinct == coinc(sel(i)));
    if down == 1
      if size(ind,2)>0
          tag(ind) = 1; % Selected candidates are tagged as valid
      else
          tag(end+1) = 1;
          coinct(end+1) = coinc(sel(i));
          runt(end+1) = run(sel(i));
      end  
    elseif down == -1
      if tag(ind) == 1
         sel2 = [sel2 sel(i)];
      end
    else
      sel2 = sel;
    end
end
[RunId ind] = sort(runt);
CoincId = coinct(ind);
  
if down==1
  sel2 = sel;
  Tag = tag(ind);
  save('CandidateTag.mat','RunId','CoincId','Tag')
  disp 'CandidateTag.mat now saved.'
else
  sel = sel2;  
end
%
disp(sprintf('%d candidates in DST, %d selected.',ncand,length(sel2)))

figure(1)
subplot(2,1,1)
hist(thetap(sel2),20)
xlim([0 90])
subplot(2,1,2)
hist(phip(sel2),40)
xlim([0 360])

figure(2)
PrepareSkyPlot(2);
polar( phip(sel2 )*DEG2RAD(1), thetap(sel2), 'go' );  

thpel = thetap(sel);
phpel = phip(sel);
radel = radius(sel);
rampel = ramp(sel);
run = run(sel);
coinc = coinc(sel);
%thsel = thetas(sel2);
%phsel = phis(sel2);
% run(sel2)
% coinc(sel2)

n115 = 0;
n137 = 0;
n128 = 0;
n112 = 0;

mult = zeros(1,length(sel));
for i=1:length(sel)
  if run(i)<3650
      continue
  end
  mult(i) = length(c.CandidateAntennas{1,sel(i)});
  antsin = c.CandidateAntennas{1,sel(i)}
  if size(find(antsin==137),2)==1
      disp 'Antenna 137 in!'
      n137 = n137+1
  end
  if size(find(antsin==115),2)==1
      disp 'Antenna 115 in!'
      n115 = n115+1
  end
  if size(find(antsin==128),2)==1
      disp 'Antenna 128 in!'
      n128 = n128+1
  end
  if size(find(antsin==112),2)==1 
      disp 'Antenna 112 in!'
      n112 = n112+1
  end  
end
n137
n115
n128
n112
pause
% 
% for i=1:length(sel)
%   CandidateTag(run(sel(i)),coinc(sel(i)))
% end

[run' coinc']

% figure(3)
% hist(c.CandidateRadius(sel),100)
