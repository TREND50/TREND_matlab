function [thpel phpel mult radel rampel run coinc val time] = ScanCandidateDST(PeriodId,down,v2)
% Adjust cuts on neighbourgs in Candidate DST
% OMH 21/06/2013
% down = 1: copy selection to CandidateTag.mat  [DANGEROUS!]
% down = 0: only displays all selected events
% down = -1: dispay selected==1 AND tag==1 events (Double selection)

SharedGlobals;
if ~exist('down')  
  down = 0;
end
if ~exist('v2')  
  v2 = 2;
end
%v2 = 1
runmax = 5914;
if v2==2
    nmax = 12;
    fend = '';
    cuts = {'L>4','Radius','Chi2','Theta<80','Barycenter','BadSignals','Pattern','Dir neighbours','Time neighbours'};
elseif v2==0
    nmax = 12;
    fend = '_check';
    cuts = {'L>4','Radius','Chi2','Theta<80','Barycenter','BadSignals','Pattern','Dir neighbours','Time neighbours'};
elseif v2==22
    nmax = 12;
    fend = '_v22';
    cuts = {'L>4','Radius','Chi2','Theta<80','Barycenter','BadSignals','Pattern','Dir neighbours','Time neighbours'};
elseif v2==3
    nmax = 13;
    fend = '_v3';
    cuts = {'L>4','Radius','Chi2','Theta<80','Silent antennas','Barycenter','BadSignals','Pattern','Dir neighbours','Time neighbours'};    
end
%     

%valdir = load('valdir_v2.txt')


%% Reads CandidateAnalysis selection results
selectionfile = sprintf('selection_Period%d.txt',PeriodId);
selectionfile = [CAND_PATH selectionfile];
a = load(selectionfile);
runid =  a(:,1);
selrun = find(runid<runmax);
% Compute number of surviving events after each cut
nsurv = zeros(size(a(:,3:nmax)));
for i = 1:length(selrun)
    % Warning: up to row 7 [theta], figures are surviving events, after
    % [n1..5] are removed events
    nsurv(i,:) = [a(i,3:7) a(i,7)-cumsum(a(i,8:nmax))];
end
% Now sum on all runs
sumsurv = sum(nsurv,1);

%badruns=[4212 4218 4271 4290 4306 4333 4352 4360]
badruns=[];

%% Load results DST
dstname = sprintf('Candidates_Period%d%s.mat',PeriodId,fend);
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
time = c.CandidateTime;
valid = zeros(1,ncand);
% ncom = 0
% for i =1:size(valdir,1)
%     if size(find(run==valdir(i,1) & coinc == valdir(i,2)),2)
%         [valdir(i,1) valdir(i,2)]
%         ncom = ncom+1
%     end
% end
% pause
minmult = 5;
n0 = 0;
n1 = 0;
n2 = 0;
n3 = 0;
n4 = 0;
n5 = 0;
n6 = 0;
n7 = 0;
n8 = 0;
fids = fopen('sandra.txt', 'a' );
    

ncand = length(find(run<runmax));
for i = 1:ncand
%     if run(i)>3601
%         break
%     end
%     if sum(ismember(badruns,run(i)))>0
%         continue
%     end
%    disp(sprintf('R%dC%d',run(i),coinc(i)))
    mult = length(c.CandidateAntennas{i});
    if mult<minmult
        n0 = n0+1;
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
    
    if thetas(i)>80 | thetas(i)<0
        %| abs(thetas(i)-thetap(i))>10 | abs(phis(i)-phip(i))>3 
        %disp 'Skip 3'
        n3 = n3+1;
        continue
    end
    
    mdirnei = dirnei{i};
%    if mdirnei(4,1)>0  % Hardest
%    if mdirnei(1,5)>0  % Softest
    if mdirnei(3,2)>0  %Best 20mn + 33% mdirnei(3,2)>0
%     if mdirnei(2,2)>0  %New best v3?
       %disp 'Skip 4.'
        n4 = n4+1;
        continue
    end
    %valdir(end+1,:)=[run(i)     coinc(i)]
    
    mnei = nei{i};
%   if mnei(4,1)>1  % Hardest 
%    if mnei(1,5)>0  %Softest
%    if mnei(2,2)>0   %New best v3?
    if mnei(1,4)>0  %Best: 30s+66% mnei(1,4)>0  
        %disp 'Skip 5.'
        n5 = n5+1;
        continue
    end
    fprintf(fids,'%10d %3d %3.2f %3.2f',time(i),run(i),thetas(i),phis(i));
    fprintf(fids,repmat(' %d',1,50),ALLDETS.*ismember(ALLDETS,ants{i}));
    fprintf(fids,'\n');
    %[time(i) run(i) thetas(i) phis(i) ALLDETS.*ismember(ALLDETS,ants{i})]    
    %pause
     if PeriodId>=9
        mants = ants{i};
         badpat = sum(ismember(mants,[135 136 137 138]));
%         badpat = sum(ismember(mants,[135 136 137 138 140 156 158]));
        if badpat>=1
            %disp 'Skip 8.'
            n8 = n8+1;
%             if sum(ismember(mants,127))>0
%               disp '127 in'
%               pause
%             end
            continue
        end
     end    
    
    %disp 'Selected.'
%     if length(ants{i}) > maxmult
%         maxmult = length(ants{i})
%         ants{i}
%         pause
%     end
    valid(i) = 1;
end
sel = find(valid==1);

val = [ncand n0 n1 n2 n6 n7 n3 n4 n5 n8 length(sel)];
nsurv2 = [ncand ncand-cumsum(val(2:10))];
% Display
display '*** Results of CandidatesAnalysis_v2 cuts ***'
for i =1:length(sumsurv)-1
    disp(sprintf('%s: %d coincs before --> %d after (%3.2f ratio)',cuts{i},sumsurv(i),sumsurv(i+1),sumsurv(i+1)./sumsurv(i)))  
end
display '*** Results of ScanCandidateDST cuts ***'
cuts2 = {'Mult','BadPulses','Chi2','Radius','Amp','Theta','Dir neighbours','Time neighbours','Flagged Antennas'};
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
time = time(sel);

for i=1:length(sel)
  if run(i)<3000  % Only consider periods 3 to 10 here
      continue
  end
  mult(i) = length(c.CandidateAntennas{1,sel(i)});
end
% selse = find(thpel>60 & phpel>210 & phpel<250)
% for i = 1:length(sel2)
%       format long
%       disp(sprintf('Coinc %d: R%dC%d, UnixSec = %d. \nAntennas: ',i,c.CandidateRun(sel2(i)),c.CandidateCoinc(sel2(i)),c.CandidateTime(sel2(i))))
%       c.CandidateAntennas{sel(i)}
% end
fclose(fids);
