function [Struct,TimeDiff,deadtime,deadtime2]=ConsecutiveCoincidenceFilter(Struct,upperlimit)


SharedGlobals;

% Load structure parameters
Evt=[Struct.Coinc.Det.Evt];
Time=[Struct.Coinc.Det.Time];
Tag=[Struct.Coinc.Det.Tag];
Ant=[Struct.Setup.Det.Name];

Reject=zeros(size(Evt,1),1);

TimeDiff=zeros(size(Evt,1),1);
EvtFlag=zeros(size(Evt,1),1);
cpt=1;

if ~exist('upperlimit')
    upperlimit=100e-3;
end;

while cpt<size(Evt,1)

    % Ref event
    Ref=cpt;
    
    % Hit antennas for ref event
    TagRef=find(Tag(Ref,:)==1);
    % Hit antennas for consecutive event
    TagCons=find(Tag(Ref+1,:)==1);
    % Minimal time on ref antenna
    TimeRef=min(Time(Ref,TagRef));
    % Minimal time on following antenna
    TimeCons=min(Time(Ref+1,TagCons));
    
    % Common antennas on both events
    ComAnt=intersect(TagRef,TagCons);
    if ~isempty(ComAnt)
         TimeDiffTotal=Time(Ref+1,ComAnt)-Time(Ref,ComAnt);
         TimeDiff(Ref)=min(TimeDiffTotal);
    else
        TimeDiffTotal=0;
        TimeDiff(Ref)=mean(Time(Ref+1,TagCons))-mean(Time(Ref,TagRef));
        %TimeDiff(Ref)=100000e9;
        EvtFlag(Ref)=1;
    end;
    clear ComAnt TagRef TagConc TimeRef TimeConc TimeDiffTotal
    
    %TimeDiff(Ref)=TimeCons-TimeRef;
    cpt=Ref+1;
    decim = floor(size(Evt,1)/10);
	if decim>0
    	if floor(cpt/decim)==cpt/decim
        	display(sprintf('Done at %2.0f percent',cpt/decim*10));
    	end
	end
end;


% To seconds
TimeDiff=TimeDiff.*5e-9;

% params rejection
freqpeak=10e-3; % peak every 10 ms
lowerlimit=1e-3; % reject all events under
%upperlimit=100e-3;
limit=3e-3;
nsigma=2;

ind=find(TimeDiff<100e-3);

Deltat=TimeDiff;

sel = find(Deltat<0.2);
figure(1);
hist(Deltat(sel)*1e3,0:0.5:200);
xlim([0 200])
xlabel('Delay (ms)')
grid on

figure(2);
[n,x] = hist(Deltat(sel)*1e3,0:2:200);
semilogy(x,n,'+k','MarkerSize',4,'LineWidth',2)
xlim([0 200])
xlabel('Delay (ms)')
grid on
pause


if ConsCoincCompleteRejection==0 
    
    indlow=find(Deltat<lowerlimit);
    Reject(indlow)=1;
    Reject(end)=0; %% last one TimeDiff=0 for sure, avoid direct rejection
    
    for i=10e-3:freqpeak:upperlimit
    
        % find maximum position for the peak
        inddt=find(Deltat<i+limit & Deltat>i-limit);
        [a,b]=hist(Deltat(inddt),min(Deltat(inddt)):0.5e-3:max(Deltat(inddt)));
        Dtmax=b(find(a==max(a)));

        if ~isempty(Dtmax)
            Dtmax=Dtmax(1);
            indmax=find(abs(Deltat-Dtmax)==min(abs(Deltat-Dtmax)));
            if length(indmax>1)
                indmax=indmax(1);
            end;
    
            %hold on
            %plot([Deltat(indmax) Deltat(indmax)],[0 10000],'r')
            % Deliminate zone around maximum
            indzone=find(Deltat<Deltat(indmax)+limit & Deltat>Deltat(indmax)-limit);
            stdzone=std(Deltat(indzone));
    
            % rejection
            indreject=find(Deltat<Deltat(indmax)+nsigma*stdzone & Deltat>Deltat(indmax)-nsigma*stdzone);
            Reject(indreject)=1;
            if indreject+1<=size(Evt,1)
                Reject(indreject+1)=1; % reject consecutive event
            end;
        end;    
        clear inddt indzone indmax Dtmax
    end;
else
    display('pouet')
    for i=1:length(TimeDiff)-1
        if TimeDiff(i)<upperlimit & EvtFlag(i)==0
            Reject(i)=1;
            Reject(i+1)=1;
        end;
    end;
end;

Struct.Coinc.Reject.ConsCoinc=Reject;

total_ant=size(Tag,2);
deadtime=0;
deadtime2=0;
% Estimation dead-time
for i=1:length(TimeDiff)-1
    TagAnt=find(Tag(i,:)==1);
    deadtime=deadtime+(length(TagAnt)/total_ant)*min(TimeDiff(i),upperlimit);
    deadtime2=deadtime2+min(TimeDiff(i),upperlimit);
end;


