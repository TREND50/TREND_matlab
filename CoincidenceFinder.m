function [Struct,RefTable,CaracDetCoinc]=CoincidenceFinder(Struct,EventTimeTable,NbIter,AnalysisType)

% Sub-structure of CoincidenceBuilder

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
sel = find(DetectorType<100);  %Anayse all detectors together
%sel=find(DetectorType==1);
Detectors = Detectors(sel);
nb_det=length(sel);
TWinFactor=1.5; % Safety factor for the time window of coincidence search (typically 1.1)

RefDet=1;
CoincNb=0;
RefTable=0;
% if scintillators exisiting, we mix the struct.setup part in one Detector
% struct
DetPosX= [RunSetup.Det.X];
DetPosY= [RunSetup.Det.Y];
DetPosZ= [RunSetup.Det.Z];

% Save tables initialisation
msize = 100000;
NbEvt=length(EventTimeTable(:,1));
TagDet=zeros(min(NbEvt,msize),length(Detectors));
StatusDet = zeros(min(NbEvt,msize),length(Detectors));
IDDet=zeros(min(NbEvt,msize),length(Detectors));
EvtDet=zeros(min(NbEvt,msize),length(Detectors));
TimeDet=zeros(min(NbEvt,msize),length(Detectors));
TriggerRate=zeros(min(NbEvt,msize),length(Detectors));
MultDet=zeros(min(NbEvt,msize),1);
MultSci=zeros(min(NbEvt,msize),1);
IDCoinc=zeros(min(NbEvt,msize),1);
UnixTime=zeros(min(NbEvt,msize),length(Detectors));

MultRejectCoinc=zeros(3,1);

display('Looking for coincidences...')
displ = 0;

while RefDet<size(EventTimeTable,1)
    
    if CoincNb>653
        displ = 0;
    end
    
    % CaracCoinc table initialisation
    % Column 1: detector ID
    % Column 2: detector type (0 : antenna, 1 : scintillator, -1 no trigger)
    % Column 3: detector event
    % Column 4: detector time trigger
    % Column 5: antenna trigger rate
    % Column 6: antenna unix time
    CaracDetCoinc = zeros(nb_det,6);
    CaracDetCoinc(:,2) = -1;
    NbDets=1;
    
    % Coincidence characteristics
    ind_refdet=find(Detectors==EventTimeTable(RefDet,2));
    PosRefDet=[DetPosX(ind_refdet) DetPosY(ind_refdet) DetPosZ(ind_refdet)];
    RefDetTime=EventTimeTable(RefDet,1);
    
    if displ
        disp 'New coinc'
        CoincNb = CoincNb
        refdet = Detectors(ind_refdet)
        EventTimeTable(RefDet,3)
    end
    
    % First detector
    CaracDetCoinc(NbDets,1)=EventTimeTable(RefDet,2);
    CaracDetCoinc(NbDets,2)=EventTimeTable(RefDet,4);
    CaracDetCoinc(NbDets,3)=EventTimeTable(RefDet,3);
    CaracDetCoinc(NbDets,4)=EventTimeTable(RefDet,1);
    CaracDetCoinc(NbDets,5)=EventTimeTable(RefDet,6);
    CaracDetCoinc(NbDets,6)=EventTimeTable(RefDet,5);
    
    for TestDet=(RefDet+1):min(RefDet+1+2*nb_det,NbEvt)
        
        TestDetTime=EventTimeTable(TestDet,1);
        ind_testdet=find(Detectors==EventTimeTable(TestDet,2));
        PosTestDet=[DetPosX(ind_testdet) DetPosY(ind_testdet) DetPosZ(ind_testdet)];
        
        if displ
            disp 'test next trig'
            Detectors(ind_testdet)
            EventTimeTable(TestDet,3)
        end
        
        % Maximal distance between the 2 studied detectors
        DistDet=norm(PosRefDet - PosTestDet);
        DistMax=DistDet*1/C0*FSAMPLING; % 1 = sin(90�) : maximal geometrical distance
        
        if displ
            TestDetTime - RefDetTime
            DistMax*TWinFactor
            pause
        end
        
        if (TestDetTime - RefDetTime)<DistMax*TWinFactor % Correction of the coincidence window with the safety factor
            % Test for double occurence of the same detector           
            if displ
                disp 'ok'
            end
            if ~isempty(find(CaracDetCoinc(:,1)==EventTimeTable(TestDet,2), 1))
                %continue
                %break; % Next coinc (switch to the next reference detector)
            else
                NbDets = NbDets+1;
                CaracDetCoinc(NbDets,1)=EventTimeTable(TestDet,2);  % DetectorId
                CaracDetCoinc(NbDets,2)=EventTimeTable(TestDet,4);  % Detector type
                CaracDetCoinc(NbDets,3)=EventTimeTable(TestDet,3);  % Evt number
                CaracDetCoinc(NbDets,4)=EventTimeTable(TestDet,1);  % Time
                CaracDetCoinc(NbDets,5)=EventTimeTable(TestDet,6);  % Trigger rate
                CaracDetCoinc(NbDets,6)=EventTimeTable(TestDet,5);  % Unix time
            end;
        else
            break; % No coinc found here, move to next one
        end;
    end;
    
    % Coincidence pre-processing
    IndScint=find(CaracDetCoinc(:,2)==1);
    NbScints = length(IndScint);  
    
    % Coincidence selection mode
    if AnalysisType==0 % radio
        CoincSelection = (NbDets>3  | NbDets==3 & NbScints>1);
    elseif AnalysisType==1 %hybrid
        CoincSelection = ((NbScints>1 & NbDets>5) | NbScints==3);
    elseif AnalysisType==2 % scint
        CoincSelection = (NbScints==3);
    end;
        
    if  CoincSelection

        CoincNb=CoincNb+1;
        
        sel = find(CaracDetCoinc(:,1)>0);
        CaracDetCoinc = CaracDetCoinc(sel,:);
        CaracDetCoinc = sortrows(CaracDetCoinc,1);
        [~,indDet]=intersect(Detectors,CaracDetCoinc(:,1));
        % 
        % Detector results saving
        MultDet(CoincNb)=NbDets;
        MultSci(CoincNb)=NbScints;
        IDCoinc(CoincNb)=msize*(NbIter-1)+CoincNb;
        IDDet(CoincNb,indDet) = CaracDetCoinc(1:NbDets,1);
        TagDet(CoincNb,indDet)=1;  % Tag detectors part of the coincidence
        StatusDet(CoincNb,indDet) = 1; % Bit0 set at alue 1. 
        EvtDet(CoincNb,indDet)=CaracDetCoinc(1:NbDets,3);
        TimeDet(CoincNb,indDet)=CaracDetCoinc(1:NbDets,4);
        UnixTime(CoincNb,indDet)=CaracDetCoinc(1:NbDets,6);
        TriggerRate(CoincNb,indDet)=CaracDetCoinc(1:NbDets,5);
        RefDet=RefDet+NbDets;
        if CoincNb>=msize
            disp(sprintf('Number of coincidences exceeds maximum allowed number (%d). Now stopping coincidence search.',msize))
            RefTable=RefDet;
            break
        end
    else
        RefDet=RefDet+1;
    end
    
    decim = floor(RunSetup.TotalEvt/10);
    if floor(RefDet/decim)==RefDet/decim
        display(sprintf('Done at %2.0f percent',RefDet/decim*10));
    end
    
    %clear CaracDetCoinc
    
end;

display(sprintf('%d coincidences found.',CoincNb))

% Mise en structure (trouver initialisation structure pour gagner du temps)
if CoincNb>0
	Struct.Coinc.Mult=MultDet(1:CoincNb);
	Struct.Coinc.MultAnt=MultDet(1:CoincNb)-MultSci(1:CoincNb);
	Struct.Coinc.MultSci=MultSci(1:CoincNb);
	Struct.Coinc.IdCoinc=IDCoinc(1:CoincNb);
	Struct.Coinc.Det.Tag=TagDet(1:CoincNb,:);
	Struct.Coinc.Det.Status=StatusDet(1:CoincNb,:);
	Struct.Coinc.Det.Id=IDDet(1:CoincNb,:);
	Struct.Coinc.Det.Evt=EvtDet(1:CoincNb,:);
	Struct.Coinc.Det.Time=TimeDet(1:CoincNb,:);
	Struct.Coinc.Det.UnixTime=UnixTime(1:CoincNb,:);
	Struct.Coinc.Det.TriggerRate=TriggerRate(1:CoincNb,:);
else
    Struct.Coinc.Mult=0;
    Struct.Coinc.MultAnt=0;
    Struct.Coinc.MultSci=0;
    Struct.Coinc.IdCoinc=0;
    Struct.Coinc.Det.Tag=0;
    Struct.Coinc.Det.Status=0;
    Struct.Coinc.Det.Id=0;
    Struct.Coinc.Det.Evt=0;
    Struct.Coinc.Det.Time=0;
    Struct.Coinc.Det.UnixTime=0;
    Struct.Coinc.Det.TriggerRate=0;
end;

Struct.Setup.TotalCoinc = CoincNb;

clear MultDet TagDet StatusDet IDDet EvtDet TimeDet;
