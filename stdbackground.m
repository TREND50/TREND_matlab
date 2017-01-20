function [Sigma,time,lst,SigmaComp,timeComp,lstComp]=stdbackground(runs,ant,comprun)
% Fetch sigma, time & lst from dst for antenna ant

SharedGlobals;
Sigma=[];
Mu=[];
MuComp=[];
Time=[];
time=[];
lst=[];
TimeComp=[];
SigmaComp=[];
timeComp=[];
lstComp=[];
nrun=[];

for i=1:length(runs)
    run = runs(i);
    disp(sprintf('stdbackground: looking for background dst for run %d...',run))
    filename=sprintf('Bgdst_%d.mat',run);
    filename=[BACKDST_PATH filename];
    ft=fopen(filename);
    if ft==-1  % Perform analysis
        disp(sprintf('stdbackground: could not find background dst for run %d, calling BackgroundDSTBuilder.',run))
        BackgroundDSTBuilder(run);
        ft=fopen(filename);
    end
    if ft~=-1  % DST was produced.
        disp(sprintf('stdbackground: found data for antenna %d in %s.',ant,filename))
        nrun=[nrun run];
        fclose(ft);
    end
end;

for i=1:length(nrun)    
    filename=sprintf('Bgdst_%d.mat',nrun(i));
    filename=[BACKDST_PATH filename];
    load(filename);
    ind=find(BgStruct.Ant==ant);
    if ~isempty(ind)
        disp(sprintf('stdbackground: uploading results of dst %d for antenna %d.',run,ant))
        Sigma=[Sigma BgStruct.Sigma(:,ind)'];
        Mu=[Mu BgStruct.Mu(:,ind)'];
        Time=[Time BgStruct.Time(:,ind)'];
        clear BgStruct;
    else
        display(sprintf('stdbackground: no data found in dst %d for antenna %d.',nrun(i),ant))
        continue
    end
    
end;

for i=1:length(Time)
    [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = unixsecs2date(Time(1,i));
    time(i)=datenum([YEAR, MONTH, DAY, HOUR, MINUTE, SECOND]);
    lst(i)=utdate2lst([YEAR, MONTH, DAY, HOUR, MINUTE, SECOND]);
end;


%% Comprun...?
if exist('comprun')
    
    for i=1:length(comprun)
        
        filename=sprintf('Bgdst_%d.mat',comprun(i));
        filename=[BACKDST_PATH filename];
        load(filename);
        ind=find(BgStruct.Ant==ant);
        if ~isempty(ind)
            SigmaComp=[SigmaComp BgStruct.Sigma(:,ind)'];
            MuComp=[MuComp BgStruct.Mu(:,ind)'];
            TimeComp=[TimeComp BgStruct.Time(:,ind)'];
            clear BgStruct;
        else
            display('Antenna not found')
            return
        end
    
    end;

    for i=1:length(TimeComp)
        [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = unixsecs2date(TimeComp(1,i));
        timeComp(i)=datenum([YEAR, MONTH, DAY, HOUR, MINUTE, SECOND]);
        lstComp(i)=utdate2lst([YEAR, MONTH, DAY, HOUR, MINUTE, SECOND]);
    
    end;
    
end;