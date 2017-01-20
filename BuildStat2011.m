function [  ] = BuildStat2011(  )
% Global stat for data run 2011
% Savec in dst data2011.mat
% OMH 10/08/2011

SharedGlobals;

scan_runs = 2550:2710;


%% load Stat2011.mat dst
if exist('Stat2011.mat')
    s = load('Stat2011.mat');
    Struct = s.Struct;
    nrun = Struct.Runs;
    totalEvts = Struct.TotalEvts;   
    date = Struct.Date;
    dur = Struct.Duration;
    detsIn = Struct.Det.Name;
    detsType = Struct.Det.isScint;
    detsEvts = Struct.Det.Evts;
    detsDur = Struct.Det.Duration;
    detsEff = Struct.Det.Eff;    
else
    nrun = [];
    totalEvts = [];
    date = [];
    dur = [];
    detsIn = [];
    detsType = [];
    detsEvts = [];
    detsDur = [];
    detsEff = [];
end

%%
for i=1:length(scan_runs)
    irun = scan_runs(i);
    if size(find(nrun==irun),2)~=0
        disp(sprintf('Run %d already in Stat2011 structure... Skipping it.',irun));
        continue
    end
    disp(sprintf('Run %d',irun));
    dstname = [DST_PATH sprintf(dst_filename,irun,1)];
    if ~exist(dstname)
        disp(sprintf('dst %s not found',dstname))
        continue
    end
    dst = load(dstname);
    nrun(end+1) = irun;
    totalEvts(end+1) = dst.Struct.Setup.TotalEvt;
    start = [dst.Struct.Setup.InfosRun.TimeStart];
    stop = [dst.Struct.Setup.InfosRun.TimeStop];
    dif = start-max(start);
    good = find(abs(dif)<24*3600);  % Keep machines within 24h of maximum TimeStart (reject machines with bad clock)
    date(end+1) = min(start(good));
    dur(end+1) = max(stop-start)/3600;  %hours

    %
    ids = [dst.Struct.Setup.Det.Name];
    detsDur(end+1,1:length(ids)) = (dst.Struct.Setup.InfosRun.TimeStop-dst.Struct.Setup.InfosRun.TimeStart)/3600;  
    detsIn(end+1,1:length(ids)) = ids;
    detsType(end+1,1:length(ids)) = [dst.Struct.Setup.Det.isScint];
    evts = [dst.Struct.Setup.Det.Evt];
    detsEvts(end+1,1:length(evts)) = evts;      
    detsEff(end+1,evts==0) = 0;
    detsEff(end,evts>0) = 1;
    detsEff(end,ids==139 | ids>=148) = detsDur(end,ids==139 | ids>=148)/dur(end);
    clear dst start stop
end

%% Load to structure
Struct.Runs = nrun;
Struct.TotalEvts = totalEvts;
Struct.Date = date;
Struct.Duration = dur;
Struct.Det.Name = detsIn;
Struct.Det.isScint = detsType;
Struct.Det.Evts = detsEvts;
Struct.Det.Duration = detsDur;
Struct.Det.Eff = detsEff;

filename = 'Stat2011.mat';
save(filename,'Struct');
clear Stat

%dst_filename = 'data2011.mat';
%if open(dst_filename)
%    dst = load(dst_filename);
%end

