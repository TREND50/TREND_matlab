function [Struct]=RunSetupBuilder(nrun,AnalysisType)
% Look for present detectors in studied run and use SetupCharacteristics to
% save importants informations as position, delays, cable...
% 14/04/10, TS
% Last modification : OMH 23/12/10

SharedGlobals;

%%%%%%%%%%%%%%% WARNING!! Change these values in case of experimental setup change %%%%%%%%%%%%%%%%%%%%%
if nrun > 444400  % NS polar
	det_search=[101:139 140 156:158 175 176];
else  %EW polar
	det_search=[101:139 140 148:158 175 176];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% SCINTILLATORS ONLY %%%%%%%%%%%%%%%%
%det_search=[139 175 176];

Detectors=[];
DetectorsFound=0;
TotalNbEvents=[];
TimeStart=[];
TimeStop=[];
DetMalfunction=[];
DetNoTrig=[];

% Identification of present antenna/scintillator in the run
for i=1:length(det_search)
    
    % Opening of antenna time file
    ft = OpenFileTime(nrun,det_search(i));
    % if exist, antenna name and number of events save
    if ft~=-1
        time=double(fread(ft,inf,'uint32'));
        fseek( ft, 0, 'eof' );
        NbEvents=ftell(ft)/(ibufft*4);
        
        if size(time,1)~=0 
            % security on acquisition bug (corrupted time file)
            unixtime(1)=time(1);
            unixtime(2)=time(1+(NbEvents-1)*ibufft);      
            if isempty(find(unixtime==0))  
                Detectors=[Detectors det_search(i)];
                DetectorsFound=DetectorsFound+1;       
                TotalNbEvents=[TotalNbEvents NbEvents]; % Total number of events for this antenna
                TimeStart=[TimeStart unixtime(1)];
                TimeStop=[TimeStop unixtime(2)];
            else
                DetMalfunction=[DetMalfunction det_search(i)];
            end
        else
            DetNoTrig=[DetNoTrig det_search(i)];
        end
        fclose(ft);
        clear unixtime NbEvents
    end;
    
end;
disp(sprintf('%d triggers on %d antennas recorded in R%d.',sum(TotalNbEvents),DetectorsFound,nrun));

% Run quality check
if (DetectorsFound==0)||(sum(TotalNbEvents)==0)
    display('RunSetupBuilder error : no access to data! Path error? ')
    return;
end;

clear RunSetup

RunSetup.Run=nrun;
RunSetup.TotalEvt=sum(TotalNbEvents);

% Setup characterisitcs (antennas, scintillators, positions, delays)

if ismember(nrun,runs2010)
    [pos_det,podN,pos_pod,cable,delay,delaycorr,sciID,logID]=SetupCharacteristics2010(Detectors);
else
    [pos_det,podN,pos_pod,cable,delay,delaycorr,sciID,logID]=SetupCharacteristics(Detectors,nrun);
end
%delay = delay-delaycorr
% RunSetup sub structre creation
cpt_det=0;

% Scintillators only analysis
if AnalysisType==2
    ind=find(sciID==1);
    Detectors=Detectors(ind);
    DetectorsFound=length(Detectors);
    TotalNbEvents=TotalNbEvents(ind);
    pos_det=pos_det(ind,:);
    podN=podN(ind);
    pos_pod=pos_pod(ind,:);
    cable=cable(ind);
    delay=delay(ind);
    delaycorr=delaycorr(ind);
    sciID=sciID(ind);
    logID=logID(ind);
    TimeStart=TimeStart(ind);
    TimeStop=TimeStop(ind);
end;

for i=1:DetectorsFound    
        cpt_det=cpt_det+1;
        det(cpt_det).Name=Detectors(i);
        det(cpt_det).Evt=TotalNbEvents(i);
        det(cpt_det).X=pos_det(i,1);
        det(cpt_det).Y=pos_det(i,2);
        det(cpt_det).Z=pos_det(i,3);
        det(cpt_det).PodNb=podN(i);
        det(cpt_det).PodX=pos_pod(i,1);
        det(cpt_det).PodY=pos_pod(i,2);
        det(cpt_det).PodZ=pos_pod(i,3);
        det(cpt_det).Cable=cable(i);
        det(cpt_det).Delay=delay(i);
        det(cpt_det).DelayCorr=delaycorr(i);
        det(cpt_det).isScint=sciID(i);     
        det(cpt_det).isLog=logID(i);
end;

% Clean TimeStart
dif = TimeStart-max(TimeStart);
good = find(abs(dif)<24*3600);  % Keep machines within 24h of maximum TimeStart (reject machines with bad clock)


Struct.Setup = RunSetup;
Struct.Setup.Det = det;
Struct.Setup.DetMalfunction=DetMalfunction;
Struct.Setup.DetNoTrig=DetNoTrig;
Struct.Setup.TotalCoinc = 0;
Struct.Setup.RunTimeStart = min(TimeStart(good));
Struct.Setup.RunTimeStop = max(TimeStop(good));
Struct.Setup.InfosRun.TimeStart = TimeStart;
Struct.Setup.InfosRun.TimeStop = TimeStop;

    
