function [  ] = AnaStat2011(  )
% Analyse dst Stat2011.mat
% with all basic global infos onacquired data 

SharedGlobals;

%% Load dst
s = load('Stat2011.mat');
Struct = s.Struct;
nrun = Struct.Runs;
totalEvts = Struct.TotalEvts;   
date = Struct.Date;
[y m d h mn s]=UnixSecs2Date(date);
dur = Struct.Duration;
detsIn = Struct.Det.Name;
detsType = Struct.Det.isScint;
detsEvts = Struct.Det.Evts;
detsDur = Struct.Det.Duration;
detsEff = Struct.Det.Eff;

nruns = length(nrun);
for i = 1:nruns
    isci = find(detsType(i,:)==1);
    iant = find(detsType(i,:)==0);
    antIn = detsIn(i,detsEvts(i,iant)>0);
    sciIn = detsIn(i,detsEvts(i,isci)>0);
    totalAnt(i) = sum(detsEvts(i,iant));
    totalSci(i) = sum(detsEvts(i,isci));
    %effAnt(i) = length(antIn)/50;
    %effSci(i) = length(sciIn)==3;
    if totalAnt(i)==0
        effAnt(i) = 0;
    else
        effAnt(i) = mean(detsEff(i,iant));
    end
    effSci(i) = mean(detsEff(i,isci))*(length(sciIn)==3);
    if ~isfinite(effSci(i))
        effSci(i) = 0;
    end
    disp(sprintf('Run %d (%02d/%02d/%d %02dh%02d) %3.1f hours\n%3.2e trigs on %d antennas (eff = %3.1f)\n%3.2e trigs on %d scints (eff=%3.1f).',nrun(i),d(i),m(i),y(i),h(i),mn(i),dur(i),totalAnt(i),length(antIn),effAnt(i)*100,totalSci(i),length(sciIn),effSci(i)*100))
end
disp ' '
%liveAnt = effAnt.*dur;
liveSci = effSci.*dur; %take out time when less than 3 scints are in...

%% Stat per month
for month = 2:4
    sel = find(m==month);
    totalAnts = sum(totalAnt(sel));
    totalScis = sum(totalSci(sel));
    totalLiveAnt = sum(dur(sel));
    totalLiveSci = sum(liveSci(sel));
    totalEffAnt = sum(effAnt(sel).*dur(sel)/sum(dur(sel)))*100;
    runsm = nrun(sel);
    disp(sprintf('Month %02d/%d (R%d-%d): \n%3.1e antenna trigs in %3.1f hours (eff = %3.1f pc, rate = %3.1fHz)\n%3.1e scint triggers in %3.1f hybrid hours (ie nscints==3), rate = %3.1f/mn',month,y(1),runsm(1),runsm(end),totalAnts,totalLiveAnt,totalEffAnt,totalAnts/totalLiveAnt/3600,totalScis,totalLiveSci,totalScis/totalLiveSci/60))
end    
figure(1)
set(1,'Name','Livetime','NumberTitle','off')

