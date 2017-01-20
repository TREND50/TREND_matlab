function [EventTimeTable,Struct]=EventTimeTableBuilder(Struct,TrigRateFilterFlag)

% Sub-structure of old CoincidenceBuilder

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
sel = find(DetectorType<100);  %Analyse all detectors together
%sel = find(DetectorType==1);
Detectors = Detectors(sel);
DetectorType=[RunSetup.Det.isScint];
Evt=[RunSetup.Det.Evt];
Evt = Evt(sel);
nb_det=length(sel);
cpt_event=0;
EventTimeTable=zeros(sum(Evt),6);

% Sort vector of all antennas events
for i=1:length(Detectors)
    
    if DetectorType(i) == 0 % antennas
        ib = ibuff;
    else
        ib = ibuffs;
    end 
    % Antenna time file reading
    ft=OpenFileTime(nrun,Detectors(i));
    time=double(fread(ft,inf,'uint32'));
    
    for j=1:Evt(i)
        cpt_event=cpt_event+1;
        ind_time=(j-1)*ibufft;
        % Time in sample since the beginning of the run
        % time(2)=nb buffer
        % time(3)=nb subbuffer
        % time(4)= position in subbuffer        
        EventTimeTable(cpt_event,1)=(time(ind_time+2)-1)*2.^28 + (time(ind_time+3)-1)*ib + time(ind_time+4) - RunSetup.Det(i).Delay;

        EventTimeTable(cpt_event,2)=Detectors(i); % Detector ID
        EventTimeTable(cpt_event,3)=j; % antenna event number
        EventTimeTable(cpt_event,4)=DetectorType(i); % if the antenna is actually a scintillator -> 1
        EventTimeTable(cpt_event,5)=time(ind_time+1); % unix seconds
       
        decim = floor((sum(Evt))/10);
        if floor(cpt_event/decim)==cpt_event/decim
            disp(sprintf('Done at %2.0f percent',cpt_event/decim*10));
        end
        
    end;
    
    clear time
    fclose(ft);
    
end;

% Remove empty signals
[EventTimeTable,Struct]=EmptySignalFilter(EventTimeTable,Struct);

% Calculations of the trigger rate for each event of each detector

DiffTime=0;
abstimerefstart=min(EventTimeTable(:,1))*(1/FSAMPLING);
abstimerefstop=max(EventTimeTable(:,1))*(1/FSAMPLING);

totaltime=abstimerefstop-abstimerefstart;
nbloop=floor(totaltime/TriggerTimeSpan);
if (totaltime-nbloop*TriggerTimeSpan)>0
    TrigRate=zeros(nbloop+1,length(Detectors));
    TrigTime=zeros(nbloop+1,1);
    DiffTime=totaltime-nbloop*TriggerTimeSpan;
else
    TrigRate=zeros(nbloop,length(Detectors));
    TrigTime=zeros(nbloop,1);
end;

time=EventTimeTable(:,1).*(1/FSAMPLING);  % in second
DetID=EventTimeTable(:,2);

for i=1:length(Detectors)
     
    display(sprintf('Trigger Rate: processing detector %d',Detectors(i)));
    inddet=find(DetID==Detectors(i));
    
    for j=1:nbloop

        ind=find(time(inddet)>abstimerefstart+(j-1)*TriggerTimeSpan & time(inddet)<abstimerefstart+j*TriggerTimeSpan);
        TrigRate(j,i)=length(ind)/TriggerTimeSpan;
        TrigTime(j)=Struct.Setup.RunTimeStart+(j-1)*TriggerTimeSpan;
        EventTimeTable(inddet(ind),6)=TrigRate(j,i);
        
        decim = floor(nbloop/10);
        if decim>0
            if floor(j/decim)==j/decim
                disp(sprintf('Done at %2.0f percent',j/decim*10));
            end
        end;
    end;
    
    if DiffTime>0
        ind=find(time>abstimerefstart+nbloop*TriggerTimeSpan & time<abstimerefstop & DetID==Detectors(i));
        TrigRate(end,i)=length(ind)/DiffTime;
        TrigTime(end)=Struct.Setup.RunTimeStart+nbloop*TriggerTimeSpan;
        EventTimeTable(ind,6)=TrigRate(end,i);
    end;
    
end;

%%%%%%%%%%%% Erase part of the event if trigger rate filtering is activated
if TrigRateFilterFlag==1
    % Dead time by antenna
    for i=1:length(Detectors)
            ind=find(TrigRate(:,i)>TrigRateLimit);
            DeadTime(i)=length(ind)*TriggerTimeSpan;
    end
    ind=find(EventTimeTable(:,6)<TrigRateLimit );  % keep scints out of this?
    display(sprintf('%d percent of the events deleted by trigger rate criterion',length(find(EventTimeTable(:,6)>=TrigRateLimit))/length(EventTimeTable(:,6))*100))
	EventTimeTable=EventTimeTable(ind,:);
end;

% Time table sort
EventTimeTable=sortrows(EventTimeTable,1);
Struct.Setup.InfosRun.TrigRate=TrigRate;
Struct.Setup.InfosRun.TrigTime=TrigTime;

if TrigRateFilterFlag==1
    Struct.Setup.InfosRun.DeadTime=DeadTime;
end;
