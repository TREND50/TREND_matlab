function [thpel phpel mult radel rampel run coinc val] = ScanCandidateDST(PeriodId,down)
% Adjust cuts on neighbourgs in Candidate DST
% OMH 21/06/2013
% down = 1: copy selection to CandidateTag.mat  [DANGEROUS!]
% down = 0: only displays all selected events
% down = -1: dispay selected==1 AND tag==1 events (Double selection)

SharedGlobals;
%     
if ~exist('down')  
  down = 0;
end

%% Reads CandidateAnalysis selection results
%selectionfile = sprintf('selection_Period%d_L4.txt',PeriodId);
selectionfile = sprintf('selection_Period%d.txt',PeriodId);
selectionfile = [CAND_PATH selectionfile];
a = load(selectionfile);
runid =  a(:,1);

% Number of surviving events after each cut
nsurv = zeros(size(a(:,3:12)));
for i = 1:size(a,1)
    % Warning: up to row 7 [theta], figures are surviving events, after
    % [n1..5] are removed events
    nsurv(i,:) = [a(i,3:7) a(i,7)-cumsum(a(i,8:12))];
%      a(i,:)
%      nsurv(i,:)
%      pause
end
% Now sum on all runs
sumsurv = sum(nsurv,1);


%badruns=[4212 4218 4271 4290 4306 4333 4352 4360]
badruns=[];

%% Load results DST
%dstname = sprintf('Candidates_Period%d_L4.mat',PeriodId);
dstname = sprintf('Candidates_Period%d.mat',PeriodId);
dstname = [CAND_PATH dstname]; 
if fopen(dstname)<0
    disp(sprintf('File %s not found. Abort.',dstname))
    thpel = [];
    phpel = [];
    mult = [];
    radel = [];
    rampel = [];
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
minmult = 5;
n1 = 0;
n2 = 0;
n3 = 0;
n4 = 0;
n5 = 0;
n6 = 0;
n7 = 0;
n8 = 0;

for i = 1:ncand
%     if sum(ismember(badruns,run(i)))>0
%         continue
%     end
    %disp(sprintf('R%dC%d',run(i),coinc(i)))
    
    mult = length(c.CandidateAntennas{i});
    if mult<minmult
        continue
    end
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
    if mdirnei(3,2)>0  %Best 20mn + 33% mdirnei(3,2)>0
       %disp 'Skip 4.'
        n4 = n4+1;
        continue
    end
    
    mnei = nei{i};
%   if mnei(4,1)>0  % Hardest 
%    if mnei(1,5)>0  %Softest
%    if mnei(3,3)>1  %Best: 2mn+50% mnei(3,3)>0   ==> ab2.ttx
    if mnei(1,4)>0  %Best: 30s+66% mnei(1,4)>0  ==> ab.txt 
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
    
%    disp 'Selected.'
%     if length(ants{i}) > maxmult
%         maxmult = length(ants{i})
%         ants{i}
%         pause
%     end
    valid(i) = 1;
end
sel = find(valid==1);
val = [ncand n1 n2 n6 n7 n3 n4 n5 n8 length(sel)];
nsurv2 = [ncand ncand-cumsum(val(2:8))];
% Display
cuts = {'L>4','Radius','Chi2','Theta<80','Barycenter','BadSignals','Pattern','Dir neighbours','Time neighbours'};
display '*** Results of CandidatesAnalysis_v2 cuts ***'
for i =1:length(sumsurv)-1
    disp(sprintf('%s: %d coincs before --> %d after (%3.2f ratio)',cuts{i},sumsurv(i),sumsurv(i+1),sumsurv(i+1)./sumsurv(i)))  
end
display '*** Results of ScanCandidateDST cuts ***'
cuts2 = {'BadPulses','Chi2','Radius','Amp','Theta','Dir neighbours','Time neighbours'};
for i = 1:length(nsurv2)-1
  disp(sprintf('%s: %d coincs before --> %d after (%3.2f ratio)',cuts2{i},nsurv2(i),nsurv2(i+1),nsurv2(i+1)./nsurv2(i)))  
end
% sumsurv'
% nsurv2'
% pause

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
disp(sprintf('***\n%d coincs, %d candidates in DSTs, %d selected (%3.1e & %3.1e ratio).',sumsurv(1),ncand,length(sel2),length(sel2)/sumsurv(1),length(sel2)/ncand))

figure(1)
subplot(2,1,1)
hist(thetap(sel2),20)
xlim([0 90])
subplot(2,1,2)
hist(phip(sel2),40)
xlim([0 360])

figure(2)
prepareSkyplot(2);
hold on
polar( phip(sel2 )*DEG2RAD(1), thetap(sel2), 'go' );  

%% Pack up returned variables
thpel = thetap(sel);
phpel = phip(sel);
radel = radius(sel);
rampel = ramp(sel);
run = run(sel);
coinc = coinc(sel);
mult = zeros(1,length(sel));
for i=1:length(sel)
  if run(i)<3650
      continue
  end
  mult(i) = length(c.CandidateAntennas{1,sel(i)});
end
  % n115 = 0;
% n137 = 0;
% n128 = 0;
% n112 = 0;
%   antsin = c.CandidateAntennas{1,sel(i)};
%   if size(find(antsin==137),2)==1
%       disp 'Antenna 137 in!'
%       n137 = n137+1;
%   end
%   if size(find(antsin==115),2)==1
%       disp 'Antenna 115 in!'
%       n115 = n115+1;
%   end
%   if size(find(antsin==128),2)==1
%       disp 'Antenna 128 in!'
%       n128 = n128+1;
%   end
%   if size(find(antsin==112),2)==1 
%       disp 'Antenna 112 in!'
%       n112 = n112+1;
%   end  
% end
