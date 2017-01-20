function []=ReducedDst(nrun)

nrun=sscanf(nrun,'%d');

SharedGlobals;

% Flag for DelayCorr fusion
DelayCorrFlag=1;


% Look for metaruns
stopflag=0;
meta=0;
AntRef=[];
while stopflag==0
    metarunname=sprintf('%04d%02d',nrun,meta);
    filename = [DST_PATH sprintf('CoincData%s_%d.mat',metarunname,1)];
    fd = fopen(filename);
    if fd~=-1
        meta=meta+1;
        fclose(fd);
        load(filename);
        AntRef=[AntRef [Struct.Setup.Det.Name]];
        clear Struct
    else
        stopflag=1;
    end;
end;
AntRef=unique(AntRef);

if meta==0
    display(sprintf('No meta run found for run %d',nrun))
    return
else
    display(sprintf('%d meta runs found for run %d',meta,nrun))
end;

% Initialisation structure variables pour Setup structure
% Setup
Run=nrun;
TotalEvt=0;
EvtAntTotal=zeros(1,length(AntRef));
DetMalfunction=[];
TotalCoinc=0;
MetaTotalCoinc=[];
RunTimeStartSave=zeros(1,meta);
RunTimeStopSave=zeros(1,meta);
SaveCoincDataTotalCoinc=zeros(1,meta);

% InfosRun
InfosSize=0;
TimeStart=zeros(1,length(AntRef));
TimeStop=zeros(1,length(AntRef));
EmptySignals=zeros(1,meta);
NbEmptySignals=0;
TrigRate=zeros(msize,length(AntRef));
TrigTime=zeros(msize,1);
DeadTime=zeros(1,length(AntRef));
GlobalCoincRateRaw=zeros(msize,1);
DetCoincRateRaw=zeros(msize,length(AntRef));
GlobalCoincRateQuickReject=zeros(msize,1);
DetCoincRateQuickReject=zeros(msize,length(AntRef));

% Coinc
Mult=[];
MultAnt=[];
MultSci=[];
IdCoinc=[];

%Coinc.Det
Tag=zeros(msize,length(AntRef));
Status=zeros(msize,length(AntRef));
Id=zeros(msize,length(AntRef));
Evt=zeros(msize,length(AntRef));
Time=zeros(msize,length(AntRef));
UnixTime=zeros(msize,length(AntRef));
TriggerRate=zeros(msize,length(AntRef));


% Finder for all referenced antennas
DetCalibSave=0;

% Preparation fusion Setup structure of meta runs
display('Fusion of Setup structures')
for i=1:meta
    
    metarunname=sprintf('%04d%02d',nrun,i-1);
    
    filename = [DST_PATH sprintf('CoincData%s_%d.mat',metarunname,1)];
    load(filename);
    display(sprintf('File %d on %d',i,meta))
    
    if Struct.Setup.TotalCoinc==0
        display('No coinc for this meta run')
        continue
    end;
    % if the metarun used all the referenced antennas, we save the Det and
    % Calib structures
    if isempty(setdiff(AntRef,[Struct.Setup.Det.Name])) & DetCalibSave==0
        Det=Struct.Setup.Det;
        Calib=Struct.Setup.InfosRun.Calib;
        DetCalibSave=1;
    end;
       
    % Start time on first file and stop time on last file
    RunTimeStartSave(i)=Struct.Setup.RunTimeStart;
    RunTimeStopSave(i)=Struct.Setup.RunTimeStop;       
    
    % Detector Malfunction
    DetMalfunction=[DetMalfunction Struct.Setup.DetMalfunction];
    
    %Empty Signals
    NbEmptySignals=NbEmptySignals+size(Struct.Setup.InfosRun.EmptySignals,1);
    EmptySignals(i)=size(Struct.Setup.InfosRun.EmptySignals,1);
    
    % In case the meta run have less antennas than the ref, we
    % determine the ind Ant index positions
    Ant=[Struct.Setup.Det.Name];
    [fake,indAnt]=intersect(AntRef,Ant);
        
    % Size InfosRun for this metarun
    InfosSizeMeta=size(Struct.Setup.InfosRun.TrigRate,1);
    % Loop on antennas
    for k=1:length(Ant)
        
        % Check for TimeStart and TimeStop
        if TimeStart(indAnt(k))==0 | TimeStart(indAnt(k))>Struct.Setup.InfosRun.TimeStart(k)
            TimeStart(indAnt(k))=Struct.Setup.InfosRun.TimeStart(k);
        end;
        if TimeStop(indAnt(k))==0 | TimeStop(indAnt(k))<Struct.Setup.InfosRun.TimeStop(k)
            TimeStop(indAnt(k))=Struct.Setup.InfosRun.TimeStop(k);
        end;
        
        % Total Evt by antenna
        EvtAntTotal(indAnt(k))=EvtAntTotal(indAnt(k))+Struct.Setup.Det(k).Evt;
        
        % InfosRun
        TrigRate(InfosSize+1:InfosSize+InfosSizeMeta,indAnt(k))=Struct.Setup.InfosRun.TrigRate(:,k);
        TrigTime(InfosSize+1:InfosSize+InfosSizeMeta)=Struct.Setup.InfosRun.TrigTime;       
        if Struct.Setup.TotalCoinc>0
            GlobalCoincRateRaw(InfosSize+1:InfosSize+InfosSizeMeta)=Struct.Setup.InfosRun.GlobalCoincRateRaw;
            DetCoincRateRaw(InfosSize+1:InfosSize+InfosSizeMeta,indAnt(k))=Struct.Setup.InfosRun.DetCoincRateRaw(:,k);
        else
            GlobalCoincRateRaw(InfosSize+1:InfosSize+InfosSizeMeta)=0;
            DetCoincRateRaw(InfosSize+1:InfosSize+InfosSizeMeta,indAnt(k))=0;
        end;   
        if isfield(Struct.Setup.InfosRun,'GlobalCoincRateQuickReject')
           GlobalCoincRateQuickReject(InfosSize+1:InfosSize+InfosSizeMeta)=Struct.Setup.InfosRun.GlobalCoincRateQuickReject;
           DetCoincRateQuickReject(InfosSize+1:InfosSize+InfosSizeMeta,indAnt(k))=Struct.Setup.InfosRun.DetCoincRateQuickReject(:,k);
        end
        if isfield(Struct.Setup.InfosRun,'DeadTime') 
           DeadTime(indAnt(k))=DeadTime(indAnt(k))+Struct.Setup.InfosRun.DeadTime(k);
        end
    end;    
    
    InfosSize=InfosSize+InfosSizeMeta;
    SaveCoincDataTotalCoinc(i)=TotalCoinc;
    TotalCoinc=TotalCoinc+Struct.Setup.TotalCoinc;
    MetaTotalCoinc=[MetaTotalCoinc Struct.Setup.TotalCoinc];
    TotalEvt=TotalEvt+Struct.Setup.TotalEvt;
  
end;

if DetCalibSave==0
    display('No meta run found with all referenced antennas, Calib and Det need recalculations')   
    [pos_det,podN,pos_pod,cable,delay,delaycorr,sciID,logID]=SetupCharacteristics(AntRef,nrun);
    for i=1:length(AntRef)    
        Det(i).Name=AntRef(i);
        Det(i).Evt=0;
        Det(i).X=pos_det(i,1);
        Det(i).Y=pos_det(i,2);
        Det(i).Z=pos_det(i,3);
        Det(i).PodNb=podN(i);
        Det(i).PodX=pos_pod(i,1);
        Det(i).PodY=pos_pod(i,2);
        Det(i).PodZ=pos_pod(i,3);
        Det(i).Cable=cable(i);
        Det(i).Delay=delay(i); 
        if DelayCorrFlag
            Det(i).DelayCorr=delaycorr(i);
        end;
        Det(i).isScint=sciID(i);     
        Det(i).isLog=logID(i);
    end;
end;

% Structure save
InfosRun.TimeStart=TimeStart;
InfosRun.TimeStop=TimeStop;
InfosRun.EmptySignals=EmptySignals;
InfosRun.MetaTotalCoinc=MetaTotalCoinc;
InfosRun.TrigRate=TrigRate(1:InfosSize,:);
InfosRun.TrigTime=TrigTime(1:InfosSize);
if DetCalibSave~=0
    InfosRun.Calib=Calib;
end;
if isfield(Struct.Setup.InfosRun,'DeadTime')
   InfosRun.DeadTime=DeadTime;
end
InfosRun.GlobalCoincRateRaw=GlobalCoincRateRaw(1:InfosSize);
InfosRun.DetCoincRateRaw=DetCoincRateRaw(1:InfosSize,:);
if isfield(Struct.Setup.InfosRun,'GlobalCoincRateQuickReject')
   InfosRun.GlobalCoincRateQuickReject=GlobalCoincRateQuickReject(1:InfosSize);
   InfosRun.DetCoincRateQuickReject=DetCoincRateQuickReject(1:InfosSize,:);
end

Setup.Run=Run;
Setup.TotalEvt=TotalEvt;
for i=1:length(AntRef)
    Det(i).Evt=EvtAntTotal(i);
end;
Setup.Det=Det;
Setup.DetMalfunction=unique(DetMalfunction);
Setup.TotalCoinc=TotalCoinc;
Setup.RunTimeStart=RunTimeStartSave;
Setup.RunTimeStop=RunTimeStopSave;
Setup.InfosRun=InfosRun;

if DetCalibSave==0;
    Struct.Setup=Setup;
    [Struct]=CalibBuilder(Struct);
    Calib=Struct.Setup.InfosRun.Calib;
    clear Struct
    Setup.InfosRun.Calib=Calib;
end;

%sum([Setup.Det.Evt])
%Setup.TotalEvt

% Reinitialisation of variabales for PreProcessData
NbBox=zeros(3*msize,length(AntRef));
TotalToT=zeros(3*msize,length(AntRef));
BoxToT=zeros(3*msize,length(AntRef));
IterId=zeros(3*msize,1);
NumCoinc=zeros(3*msize,1);

display('Recup TimeDiff informations...')
% Recup Time diff sur fichier coinc
% Time diff variables
IterIdDt=zeros(TotalCoinc,1);
NumCoincDt=zeros(TotalCoinc,1);
TimeDiffTotal=zeros(TotalCoinc,1);
NbCoincFinal=0;

for i=1:meta
    
    metarunname=sprintf('%04d%02d',nrun,i-1);
    display(sprintf('File %d on %d',i,meta))
    %Looking for substructure
    stopflag=0;
    nbiter=2;
    while stopflag==0
        filename = [DST_PATH sprintf('CoincData%s_%d.mat',metarunname,nbiter)];
        fd=fopen(filename);
        if fd~=-1
            nbiter=nbiter+1;
            fclose(fd);
        else
            stopflag=1;
        end;
    end;
    nbiter=nbiter-1;
    
    % Loop on substructure   
    for j=1:nbiter
        
        filename = [DST_PATH sprintf('CoincData%s_%d.mat',metarunname,j)];
        load(filename);
        if Struct.Setup.TotalCoinc==0
            display(sprintf('No coinc for meta run %s iteration %d',metarunname,j))
            continue
        end;
        nbcoinc=Struct.Setup.TotalCoinc;
        [~,TimeDiff]=ConsecutiveCoincidenceFilter(Struct);
        
        for k=1:nbcoinc
            
            NbCoincFinal=NbCoincFinal+1;
            IterIdDt(NbCoincFinal)=i;
            NumCoincDt(NbCoincFinal)=Struct.Coinc.IdCoinc(k);
            TimeDiffTotal(NbCoincFinal)=TimeDiff(k);
            
        end;
        clear TimeDiff
        
    end;
end;
           
% Fusion PreProcessDataFile
display('Informations of PreProcessData structures')

% Check for msize limit
NbCoincFinal=0;
IterFinal=1;
EvtAntTotal=zeros(1,length(AntRef));
TotalCoinc=0;

for i=1:meta
    
    metarunname=sprintf('%04d%02d',nrun,i-1);
    display(sprintf('File %d on %d',i,meta))
    %Looking for substructure
    stopflag=0;
    nbiter=2;
    while stopflag==0
        filename = [DST_PATH sprintf('PreProcessData%s_%d.mat',metarunname,nbiter)];
        fd=fopen(filename);
        if fd~=-1
            nbiter=nbiter+1;
            fclose(fd);
        else
            stopflag=1;
        end;
    end;
    nbiter=nbiter-1;
    
    % Loop on substructure   
    CoincMetaIter=0;
    for j=1:nbiter
        
        if j>1
	    clear StructT
        end;
        filename = [DST_PATH sprintf('PreProcessData%s_%d.mat',metarunname,j)];
        load(filename);
        if Struct.Setup.TotalCoinc==0
            display(sprintf('No coinc for meta run %s iteration %d',metarunname,j))
            continue
        end;
        StructT=Struct; % temporary version of Struct
        clear Struct
        
        % Identification of antenas in metarun
        Ant=[StructT.Setup.Det.Name];
        [fake,indAnt]=intersect(AntRef,Ant);
       
        CoincMeta=StructT.Setup.TotalCoinc;
        TagPreProcess=StructT.Coinc.Det.Tag;
        IdCoincPreProcess=StructT.Coinc.IdCoinc;
        CaracEvt=StructT.Coinc.Reject.CaracEvt;
        %NbBoxStruct=[CaracEvt{:,:,2}];
        %TotalToTStruct=[CaracEvt{:,:,3}];
        %bstartStruct=[CaracEvt{:,:,4}];
        %bstopStruct=[CaracEvt{:,:,5}];
        %blengthStruct=[CaracEvt{:,:,6}];

       tic 
        % Loop on CoincMeta
        for k=1:CoincMeta
            
            NbCoincFinal=NbCoincFinal+1;
            IterId(NbCoincFinal)=i;
            NumCoinc(NbCoincFinal)=IdCoincPreProcess(k);
            
            ind=find(TagPreProcess(k,:)==1);
		    
            NbBox(NbCoincFinal,indAnt(ind))=[CaracEvt{k,ind,2}];
            TotalToT(NbCoincFinal,indAnt(ind))=[CaracEvt{k,ind,3}];
            bstart=[CaracEvt{k,ind,4}];
            bstop=[CaracEvt{k,ind,5}];
            blength=[CaracEvt{k,ind,6}];
            start=1;
            for m=1:length(ind)
                antstart=bstart([start:start+NbBox(NbCoincFinal,indAnt(ind(m)))-1]);
                antstop=bstop([start:start+NbBox(NbCoincFinal,indAnt(ind(m)))-1]);
                indcenter=find(antstart<ibuff/2 & antstop>ibuff/2);
                if ~isempty(indcenter)
                    BoxToT(NbCoincFinal,indAnt(ind(m)))=blength(start+indcenter-1);
                end
                start=start+NbBox(NbCoincFinal,indAnt(ind(m)));
            end;
            clear bstart bstop indcenter blength start antstart antstop
		                   
        end;
        toc
        clear TagPreProcess IdCoinc CaracEvt      
    end;
end;

clear StructT

IterId=IterId(1:NbCoincFinal);
NumCoinc=NumCoinc(1:NbCoincFinal);                
NbBox=NbBox(1:NbCoincFinal,:);
TotalToT=TotalToT(1:NbCoincFinal,:);
BoxToT=BoxToT(1:NbCoincFinal,:);

% Initialisation variables for dst file fusion
% Coinc
Mult=[];
MultAnt=[];
IdCoinc=[];
IsShower=[];

% Coinc.Det
Status=-1*ones(msize,length(AntRef));
Evt=zeros(msize,length(AntRef));
UnixTime=zeros(msize,length(AntRef));
TriggerRate=zeros(msize,length(AntRef));
Sigma=zeros(msize,length(AntRef));
MinRaw=zeros(msize,length(AntRef));
MaxRaw=zeros(msize,length(AntRef));
AmpMax=zeros(msize,length(AntRef));
TrigTime=zeros(msize,length(AntRef));
Gain=zeros(msize,length(AntRef));
TrigCor=zeros(msize,length(AntRef));
CoefCor=zeros(msize,length(AntRef));
GainPSD=cell(1,msize);
CalibratedAmp1=zeros(msize,length(AntRef));
CalibratedAmp2=zeros(msize,length(AntRef));
Data3Scints=cell(msize,length(AntRef));

% Coinc.PlanRecons

% Radio
FlagR=zeros(1,msize);
ThetaR=zeros(1,msize);
dThetaR=zeros(1,msize);
PhiR=zeros(1,msize);
dPhiR=zeros(1,msize);
Chi2R=zeros(1,msize);
SignifR=zeros(1,msize);
Chi2DelayR=zeros(1,msize);
SlopeDelayR=zeros(1,msize);

% Hybrid
FlagH=zeros(1,msize);
ThetaH=zeros(1,msize);
dThetaH=zeros(1,msize);
PhiH=zeros(1,msize);
dPhiH=zeros(1,msize);
Chi2H=zeros(1,msize);
SignifH=zeros(1,msize);
Chi2DelayH=zeros(1,msize);
SlopeDelayH=zeros(1,msize);

% Sci
FlagS=zeros(1,msize);
ThetaS=zeros(1,msize);
dThetaS=zeros(1,msize);
PhiS=zeros(1,msize);
dPhiS=zeros(1,msize);

% Coinc.SphRecons
FlagSp=zeros(1,msize);
Rho=zeros(1,msize);
ThetaSp=zeros(1,msize);
PhiSp=zeros(1,msize);
X0=zeros(1,msize);
Y0=zeros(1,msize);
Z0=zeros(1,msize);
DistSource=zeros(msize,length(AntRef));
Chi2DelaySp=zeros(1,msize);
SlopeDelaySp=zeros(1,msize);

if DelayCorrFlag
    % Coinc.DelayCorr.PlanRecons
    FlagCorr=zeros(1,msize);
    ThetaCorr=zeros(1,msize);
    dThetaCorr=zeros(1,msize);
    PhiCorr=zeros(1,msize);
    dPhiCorr=zeros(1,msize);
    Chi2Corr=zeros(1,msize);
    SignifCorr=zeros(1,msize);
    Chi2DelayCorr=zeros(1,msize);
    SlopeDelayCorr=zeros(1,msize);

    % Coinc.DelayCorr.SphRecons
    FlagSpCorr=zeros(1,msize);
    RhoCorr=zeros(1,msize);
    ThetaSpCorr=zeros(1,msize);
    PhiSpCorr=zeros(1,msize);
    X0Corr=zeros(1,msize);
    Y0Corr=zeros(1,msize);
    Z0Corr=zeros(1,msize);
    DistSourceCorr=zeros(msize,length(AntRef));
    Chi2DelaySpCorr=zeros(1,msize);
    SlopeDelaySpCorr=zeros(1,msize);
end;

% Coinc.ShowerRecons
ShowerSignals=cell(msize,length(AntRef));
XCore=zeros(msize,1);
YCore=zeros(msize,1);
ZCore=zeros(msize,1);
AxisAmp=zeros(msize,1);
Lambda=zeros(msize,1);

% Coinc.InfosSignal
NbBoxSave=zeros(msize,length(AntRef));
TotalToTSave=zeros(msize,length(AntRef));
BoxToTSave=zeros(msize,length(AntRef));
TimeDiffSave=zeros(1,msize);


% Fusion DST File
display('Fusion of DST structures')


% Check for msize limit
NbCoincFinal=0;
IterFinal=1;
EvtAntTotal=zeros(1,length(AntRef));
TotalCoinc=0;


for i=1:meta
    
    metarunname=sprintf('%04d%02d',nrun,i-1);
    display(sprintf('File %d on %d',i,meta))
    %Looking for substructure
    stopflag=0;
    nbiter=2;
    while stopflag==0
        filename = [DST_PATH sprintf('dst%s_%d.mat',metarunname,nbiter)];
        fd=fopen(filename);
        if fd~=-1
            nbiter=nbiter+1;
            fclose(fd);
        else
            stopflag=1;
        end;
    end;
    nbiter=nbiter-1;
    
    % Loop on substructure 
    CoincMetaIter=0;
    for j=1:nbiter
        
	if j>1
	    clear StructT
	end;
        filename = [DST_PATH sprintf('dst%s_%d.mat',metarunname,j)];
        load(filename);
        
        % Identification of antenas in metarun
        Ant=[StructT.Setup.Det.Name];
        [fake,indAnt]=intersect(AntRef,Ant);
        
        if Struct.Setup.TotalCoinc==0
            display(sprintf('No coinc for meta run %s iteration %d',metarunname,j))
            continue
        end;
        [Struct]=Dist2Source(Struct);
        StructT=Struct; % temporary version of Struct
        clear Struct
        
        CoincMeta=StructT.Setup.TotalCoinc;
        
        % Loop on CoincMeta
        for k=1:CoincMeta
            
            NbCoincFinal=NbCoincFinal+1;
            
            if NbCoincFinal<=msize
                
                % Recup time diff
                indcoincdt=find(IterIdDt==i & NumCoincDt==StructT.Coinc.IdCoinc(k));
                if ~isempty(indcoincdt)
                    TimeDiffSave(NbCoincFinal)=TimeDiffTotal(indcoincdt);
                end;
                
                % Recup signal informations
                indcoinc=find(IterId==i & NumCoinc==StructT.Coinc.IdCoinc(k));
                if ~isempty(indcoinc)
                    NbBoxSave(NbCoincFinal,:)=NbBox(indcoinc,:);
                    TotalToTSave(NbCoincFinal,:)=TotalToT(indcoinc,:);
                    BoxToTSave(NbCoincFinal,:)=BoxToT(indcoinc,:);
                end;
                
                Mult(NbCoincFinal)=StructT.Coinc.Mult(k);
                MultAnt(NbCoincFinal)=StructT.Coinc.MultAnt(k);
                IdCoinc(NbCoincFinal)=StructT.Coinc.IdCoinc(k)+SaveCoincDataTotalCoinc(i); %shift of IdCoinc depending on TotalCoinc
                IsShower(NbCoincFinal)=StructT.Coinc.IsShower(k);
                
                Status(NbCoincFinal,indAnt)=StructT.Coinc.Det.Status(k,:);
                UnixTime(NbCoincFinal,indAnt)=StructT.Coinc.Det.UnixTime(k,:);
                TriggerRate(NbCoincFinal,indAnt)=StructT.Coinc.Det.TriggerRate(k,:);
                Sigma(NbCoincFinal,indAnt)=StructT.Coinc.Det.Sigma(k,:);
                MinRaw(NbCoincFinal,indAnt)=StructT.Coinc.Det.MinRaw(k,:);
                MaxRaw(NbCoincFinal,indAnt)=StructT.Coinc.Det.MaxRaw(k,:);
                AmpMax(NbCoincFinal,indAnt)=StructT.Coinc.Det.AmpMax(k,:);
                TrigTime(NbCoincFinal,indAnt)=StructT.Coinc.Det.TrigTime(k,:);
                Gain(NbCoincFinal,indAnt)=StructT.Coinc.Det.Gain(k,:);
                TrigCor(NbCoincFinal,indAnt)=StructT.Coinc.Det.TrigCor(k,:);
                CoefCor(NbCoincFinal,indAnt)=StructT.Coinc.Det.CoefCor(k,:);
                GainPSD(NbCoincFinal)=StructT.Coinc.Det.GainPSD(k);
                CalibratedAmp1(NbCoincFinal,indAnt)=StructT.Coinc.Det.CalibratedAmp1(k,:);
                CalibratedAmp2(NbCoincFinal,indAnt)=StructT.Coinc.Det.CalibratedAmp2(k,:);
                for m=1:length(Ant)
                    if StructT.Coinc.Det.Tag(k,m)==1
                        Evt(NbCoincFinal,indAnt(m))=StructT.Coinc.Det.Evt(k,m)+EvtAntTotal(indAnt(m));
                        Data3Scints(NbCoincFinal,indAnt(m))=StructT.Coinc.Det.Data3Scints(k,m);
                    end;
                end;
                
                %Radio
                FlagR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Flag(k);
                ThetaR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Theta(k);
                dThetaR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.dTheta(k);
                PhiR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Phi(k);
                dPhiR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.dPhi(k);
                Chi2R(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Chi2(k);
                SignifR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Signif(k);
                Chi2DelayR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Chi2Delay(k);
                SlopeDelayR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.SlopeDelay(k);
                
                % Hybrid
                FlagH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Flag(k);
                ThetaH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Theta(k);
                dThetaH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.dTheta(k);
                PhiH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Phi(k);
                dPhiH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.dPhi(k);
                Chi2H(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Chi2(k);
                SignifH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Signif(k);
                Chi2DelayH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Chi2Delay(k);
                SlopeDelayH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.SlopeDelay(k);
                
                % Sci
                FlagS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Flag(k);
                ThetaS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Theta(k);
                dThetaS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.dTheta(k);
                PhiS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Phi(k);
                dPhiS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.dPhi(k);

                % Coinc.SphRecons
                FlagSp(NbCoincFinal)=StructT.Coinc.SphRecons.Flag(k);
                Rho(NbCoincFinal)=StructT.Coinc.SphRecons.Rho(k);
                ThetaSp(NbCoincFinal)=StructT.Coinc.SphRecons.Theta(k);
                PhiSp(NbCoincFinal)=StructT.Coinc.SphRecons.Phi(k);
                X0(NbCoincFinal)=StructT.Coinc.SphRecons.X0(k);
                Y0(NbCoincFinal)=StructT.Coinc.SphRecons.Y0(k);
                Z0(NbCoincFinal)=StructT.Coinc.SphRecons.Z0(k);
                DistSource(NbCoincFinal,indAnt)=StructT.Coinc.SphRecons.DistSource(k,:);
                Chi2DelaySp(NbCoincFinal)=StructT.Coinc.SphRecons.Chi2Delay(k);
                SlopeDelaySp(NbCoincFinal)=StructT.Coinc.SphRecons.SlopeDelay(k);
                
                if DelayCorrFlag
                    % Coinc.DelayCorr.PlanRecons
                    FlagCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Flag(k);
                    ThetaCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Theta(k);
                    dThetaCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.dTheta(k);
                    PhiCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Phi(k);
                    dPhiCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.dPhi(k);
                    Chi2Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Chi2(k);
                    SignifCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Signif(k);
                    Chi2DelayCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Chi2Delay(k);
                    SlopeDelayCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.SlopeDelay(k);
                
                    % Coinc.DelayCorr.SphRecons
                    FlagSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Flag(k);
                    RhoCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Rho(k);
                    ThetaSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Theta(k);
                    PhiSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Phi(k);
                    X0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.X0(k);
                    Y0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Y0(k);
                    Z0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Z0(k);
                    DistSourceCorr(NbCoincFinal,indAnt)=StructT.Coinc.DelayCorrRecons.SphRecons.DistSource(k,:);
                    Chi2DelaySpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Chi2Delay(k);
                    SlopeDelaySpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.SlopeDelay(k);
                end
                
                % Coinc.ShowerRecons
                for m=1:length(Ant)
                    ShowerSignals(NbCoincFinal,indAnt(m))=StructT.Coinc.ShowerRecons.ShowerSignals(k,m);
                end;
                if isfield(StructT.Coinc.ShowerRecons,'XCore')
                    XCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.XCore(k);
                    YCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.YCore(k);
                    ZCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.ZCore(k);
                    AxisAmp(NbCoincFinal)=StructT.Coinc.ShowerRecons.AxisAmp(k);
                    Lambda(NbCoincFinal)=StructT.Coinc.ShowerRecons.Lambda(k);
                end;
                
            else   % we have to save CoincData and move to next iteration
            
                Struct.Setup=Setup;
                Struct.Setup.TotalCoinc=NbCoincFinal-1; % need to find a way to recover TotalCoinc in that case
                Struct.Coinc.Mult=Mult';
                Struct.Coinc.MultAnt=MultAnt';
                Struct.Coinc.IdCoinc=IdCoinc';
                Struct.Coinc.IsShower=IsShower';
                
                Struct.Coinc.Det.Status=Status;
                Struct.Coinc.Det.Evt=Evt;
                Struct.Coinc.Det.UnixTime=UnixTime;
                Struct.Coinc.Det.TriggerRate=TriggerRate;
                Struct.Coinc.Det.Sigma=Sigma;
                Struct.Coinc.Det.MinRaw=MinRaw;
                Struct.Coinc.Det.MaxRaw=MaxRaw;
                Struct.Coinc.Det.AmpMax=AmpMax;
                Struct.Coinc.Det.TrigTime=TrigTime;
                Struct.Coinc.Det.Gain=Gain;
                Struct.Coinc.Det.TrigCor=TrigCor;
                Struct.Coinc.Det.CoefCor=CoefCor;
                Struct.Coinc.Det.GainPSD=GainPSD;
                Struct.Coinc.Det.CalibratedAmp1=CalibratedAmp1;
                Struct.Coinc.Det.CalibratedAmp2=CalibratedAmp2;
                Struct.Coinc.Det.Data3Scints=Data3Scints;
                
                Struct.Coinc.PlanRecons.Radio.Flag=FlagR;
                Struct.Coinc.PlanRecons.Radio.Theta=ThetaR;
                Struct.Coinc.PlanRecons.Radio.dTheta=dThetaR;
                Struct.Coinc.PlanRecons.Radio.Phi=PhiR;
                Struct.Coinc.PlanRecons.Radio.dPhi=dPhiR;
                Struct.Coinc.PlanRecons.Radio.Chi2=Chi2R;
                Struct.Coinc.PlanRecons.Radio.Signif=SignifR;
                Struct.Coinc.PlanRecons.Radio.Chi2Delay=Chi2DelayR;
                Struct.Coinc.PlanRecons.Radio.SlopeDelay=SlopeDelayR;
                
                Struct.Coinc.PlanRecons.Hybrid.Flag=FlagH;
                Struct.Coinc.PlanRecons.Hybrid.Theta=ThetaH;
                Struct.Coinc.PlanRecons.Hybrid.dTheta=dThetaH;
                Struct.Coinc.PlanRecons.Hybrid.Phi=PhiH;
                Struct.Coinc.PlanRecons.Hybrid.dPhi=dPhiH;
                Struct.Coinc.PlanRecons.Hybrid.Chi2=Chi2H;
                Struct.Coinc.PlanRecons.Hybrid.Signif=SignifH;
                Struct.Coinc.PlanRecons.Hybrid.Chi2Delay=Chi2DelayH;
                Struct.Coinc.PlanRecons.Hybrid.SlopeDelay=SlopeDelayH;
                
                Struct.Coinc.PlanRecons.Sci.Flag=FlagS;
                Struct.Coinc.PlanRecons.Sci.Theta=ThetaS;
                Struct.Coinc.PlanRecons.Sci.dTheta=dThetaS;
                Struct.Coinc.PlanRecons.Sci.Phi=PhiS;
                Struct.Coinc.PlanRecons.Sci.dPhi=dPhiS;
                
                Struct.Coinc.SphRecons.Flag=FlagSp;
                Struct.Coinc.SphRecons.Rho=Rho;
                Struct.Coinc.SphRecons.Theta=ThetaSp;
                Struct.Coinc.SphRecons.Phi=PhiSp;
                Struct.Coinc.SphRecons.X0=X0;
                Struct.Coinc.SphRecons.Y0=Y0;
                Struct.Coinc.SphRecons.Z0=Z0;
                Struct.Coinc.SphRecons.DistSource=DistSource;
                Struct.Coinc.SphRecons.Chi2Delay=Chi2DelaySp;
                Struct.Coinc.SphRecons.SlopeDelay=SlopeDelaySp;
                
                if DelayCorrFlag
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Flag=FlagCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Theta=ThetaCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.dTheta=dThetaCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Phi=PhiCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.dPhi=dPhiCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Chi2=Chi2Corr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Signif=SignifCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.Chi2Delay=Chi2DelayCorr;
                    Struct.Coinc.DelayCorrRecons.PlanRecons.SlopeDelay=SlopeDelayCorr;
                
                    Struct.Coinc.DelayCorrRecons.SphRecons.Flag=FlagSpCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Rho=RhoCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Theta=ThetaSpCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Phi=PhiSpCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.X0=X0Corr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Y0=Y0Corr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Z0=Z0Corr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.DistSource=DistSourceCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.Chi2Delay=Chi2DelaySpCorr;
                    Struct.Coinc.DelayCorrRecons.SphRecons.SlopeDelay=SlopeDelaySpCorr;
                end
                
                Struct.Coinc.ShowerRecons.ShowerSignals=ShowerSignals;
                if isfield(StructT.Coinc.ShowerRecons,'XCore')
                    Struct.Coinc.ShowerRecons.XCore=XCore;
                    Struct.Coinc.ShowerRecons.YCore=YCore;
                    Struct.Coinc.ShowerRecons.ZCore=ZCore;
                	Struct.Coinc.ShowerRecons.AxisAmp=AxisAmp;
                    Struct.Coinc.ShowerRecons.Lambda=Lambda;
                end;
                
                Struct.Coinc.InfosSignal.NbBox=NbBoxSave;
                Struct.Coinc.InfosSignal.TotalToT=TotalToTSave;
                Struct.Coinc.InfosSignal.BoxToT=BoxToTSave;
                Struct.Coinc.InfosSignal.TimeDiff=TimeDiffSave;
                
                filename=[DST_PATH sprintf(dst_filename,nrun,IterFinal)];
                save(filename,'Struct');
                clear Struct
                IterFinal=IterFinal+1;
                NbCoincFinal=1;
                
                % Reinitialisation CoincData variables          
                Mult=[];
                MultAnt=[];
                IdCoinc=[];
                IsShower=[];
                Status=-1*ones(msize,length(AntRef));
                Evt=zeros(msize,length(AntRef));
                UnixTime=zeros(msize,length(AntRef));
                TriggerRate=zeros(msize,length(AntRef));
                Sigma=zeros(msize,length(AntRef));
                MinRaw=zeros(msize,length(AntRef));
                MaxRaw=zeros(msize,length(AntRef));
                AmpMax=zeros(msize,length(AntRef));
                TrigTime=zeros(msize,length(AntRef));
                Gain=zeros(msize,length(AntRef));
                TrigCor=zeros(msize,length(AntRef));
                CoefCor=zeros(msize,length(AntRef));
                GainPSD=cell(1,msize);
                CalibratedAmp1=zeros(msize,length(AntRef));
                CalibratedAmp2=zeros(msize,length(AntRef));
                Data3Scints=cell(msize,length(AntRef));

                % Coinc.PlanRecons
                % Radio
                FlagR=zeros(1,msize);
                ThetaR=zeros(1,msize);
                dThetaR=zeros(1,msize);
                PhiR=zeros(1,msize);
                dPhiR=zeros(1,msize);
                Chi2R=zeros(1,msize);
                SignifR=zeros(1,msize);
                Chi2DelayR=zeros(1,msize);
                SlopeDelayR=zeros(1,msize);

                % Hybrid
                FlagH=zeros(1,msize);
                ThetaH=zeros(1,msize);
                dThetaH=zeros(1,msize);
                PhiH=zeros(1,msize);
                dPhiH=zeros(1,msize);
                Chi2H=zeros(1,msize);
                SignifH=zeros(1,msize);
                Chi2DelayH=zeros(1,msize);
                SlopeDelayH=zeros(1,msize);

                % Sci
                FlagS=zeros(1,msize);
                ThetaS=zeros(1,msize);
                dThetaS=zeros(1,msize);
                PhiS=zeros(1,msize);
                dPhiS=zeros(1,msize);

                % Coinc.SphRecons
                FlagSp=zeros(1,msize);
                Rho=zeros(1,msize);
                ThetaSp=zeros(1,msize);
                PhiSp=zeros(1,msize);
                X0=zeros(1,msize);
                Y0=zeros(1,msize);
                Z0=zeros(1,msize);
                DistSource=zeros(msize,length(AntRef));
                Chi2DelaySp=zeros(1,msize);
                SlopeDelaySp=zeros(1,msize);
                
                if DelayCorrFlag
                    % Coinc.DelayCorr.PlanRecons
                    FlagCorr=zeros(1,msize);
                    ThetaCorr=zeros(1,msize);
                    dThetaCorr=zeros(1,msize);
                    PhiCorr=zeros(1,msize);
                    dPhiCorr=zeros(1,msize);
                    Chi2Corr=zeros(1,msize);
                    SignifCorr=zeros(1,msize);
                    Chi2DelayCorr=zeros(1,msize);
                    SlopeDelayCorr=zeros(1,msize);

                    % Coinc.DelayCorr.SphRecons
                    FlagSpCorr=zeros(1,msize);
                    RhoCorr=zeros(1,msize);
                    ThetaSpCorr=zeros(1,msize);
                    PhiSpCorr=zeros(1,msize);
                    X0Corr=zeros(1,msize);
                    Y0Corr=zeros(1,msize);
                    Z0Corr=zeros(1,msize);
                    DistSourceCorr=zeros(msize,length(AntRef));
                    Chi2DelaySpCorr=zeros(1,msize);
                    SlopeDelaySpCorr=zeros(1,msize);
                end;

                % Coinc.ShowerRecons
                ShowerSignals=cell(msize,length(AntRef));
                XCore=zeros(msize,1);
                Ycore=zeros(msize,1);
                ZCore=zeros(msize,1);
                AxisAmp=zeros(msize,1);
                Lambda=zeros(msize,1);
                
                % Coinc.InfosSignal
                NbBoxSave=zeros(msize,length(AntRef));
                TotalToTSave=zeros(msize,length(AntRef));
                BoxToTSave=zeros(msize,length(AntRef));
                TimeDiffSave=zeros(1,msize);
                
                % Save the missing coinc
                Mult(NbCoincFinal)=StructT.Coinc.Mult(k);
                MultAnt(NbCoincFinal)=StructT.Coinc.MultAnt(k);
                IdCoinc(NbCoincFinal)=StructT.Coinc.IdCoinc(k)+TotalCoinc; %shift of IdCoinc depending on TotalCoinc
                IsShower(NbCoincFinal)=StructT.Coinc.IsShower(k);
                
                % Recup time diff
                indcoincdt=find(IterIdDt==i & NumCoincDt==StructT.Coinc.IdCoinc(k));
                if ~isempty(indcoincdt)
                    TimeDiffSave(NbCoincFinal)=TimeDiffTotal(indcoincdt);
                end;
                
                % Recup signal informations
                indcoinc=find(IterId==j & NumCoinc==StructT.Coinc.IdCoinc(NbCoincFinal));
                if isempty(indcoinc)
                    NbBoxSave(NbCoincFinal,:)=NbBox(indcoinc,:);
                    TotalToTSave(NbCoincFinal,:)=TotalToT(indcoinc,:);
                    BoxToTSave(NbCoincFinal,:)=BoxToT(indcoinc,:);
                end;
                
                Status(NbCoincFinal,indAnt)=StructT.Coinc.Det.Status(k,:);
                UnixTime(NbCoincFinal,indAnt)=StructT.Coinc.Det.UnixTime(k,:);
                TriggerRate(NbCoincFinal,indAnt)=StructT.Coinc.Det.TriggerRate(k,:);
                Sigma(NbCoincFinal,indAnt)=StructT.Coinc.Det.Sigma(k,:);
                MinRaw(NbCoincFinal,indAnt)=StructT.Coinc.Det.MinRaw(k,:);
                MaxRaw(NbCoincFinal,indAnt)=StructT.Coinc.Det.MaxRaw(k,:);
                AmpMax(NbCoincFinal,indAnt)=StructT.Coinc.Det.AmpMax(k,:);
                TrigTime(NbCoincFinal,indAnt)=StructT.Coinc.Det.TrigTime(k,:);
                Gain(NbCoincFinal,indAnt)=StructT.Coinc.Det.Gain(k,:);
                TrigCor(NbCoincFinal,indAnt)=StructT.Coinc.Det.TrigCor(k,:);
                CoefCor(NbCoincFinal,indAnt)=StructT.Coinc.Det.CoefCor(k,:);
                GainPSD(NbCoincFinal)=StructT.Coinc.Det.GainPSD(k);
                CalibratedAmp1(NbCoincFinal,indAnt)=StructT.Coinc.Det.CalibratedAmp1(k,:);
                CalibratedAmp2(NbCoincFinal,indAnt)=StructT.Coinc.Det.CalibratedAmp2(k,:);
                for m=1:length(Ant)
                    if StructT.Coinc.Det.Tag(k,m)==1
                        Evt(NbCoincFinal,indAnt(m))=StructT.Coinc.Det.Evt(k,m)+EvtAntTotal(indAnt(m));
                        Data3Scints(NbCoincFinal,indAnt(m))=StructT.Coinc.Det.Data3Scints(k,m);
                    end;
                end;
                
                %Radio
                FlagR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Flag(k);
                ThetaR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Theta(k);
                dThetaR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.dTheta(k);
                PhiR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Phi(k);
                dPhiR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.dPhi(k);
                Chi2R(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Chi2(k);
                SignifR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Signif(k);
                Chi2DelayR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.Chi2Delay(k);
                SlopeDelayR(NbCoincFinal)=StructT.Coinc.PlanRecons.Radio.SlopeDelay(k);

                % Hybrid
                FlagH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Flag(k);
                ThetaH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Theta(k);
                dThetaH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.dTheta(k);
                PhiH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Phi(k);
                dPhiH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.dPhi(k);
                Chi2H(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Chi2(k);
                SignifH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Signif(k);
                Chi2DelayH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.Chi2Delay(k);
                SlopeDelayH(NbCoincFinal)=StructT.Coinc.PlanRecons.Hybrid.SlopeDelay(k);
                
                % Sci
                FlagS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Flag(k);
                ThetaS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Theta(k);
                dThetaS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.dTheta(k);
                PhiS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.Phi(k);
                dPhiS(NbCoincFinal)=StructT.Coinc.PlanRecons.Sci.dPhi(k);

                % Coinc.SphRecons
                FlagSp(NbCoincFinal)=StructT.Coinc.SphRecons.Flag(k);
                Rho(NbCoincFinal)=StructT.Coinc.SphRecons.Rho(k);
                ThetaSp(NbCoincFinal)=StructT.Coinc.SphRecons.Theta(k);
                PhiSp(NbCoincFinal)=StructT.Coinc.SphRecons.Phi(k);
                X0(NbCoincFinal)=StructT.Coinc.SphRecons.X0(k);
                Y0(NbCoincFinal)=StructT.Coinc.SphRecons.Y0(k);
                Z0(NbCoincFinal)=StructT.Coinc.SphRecons.Z0(k);
                DistSource(NbCoincFinal,indAnt)=StructT.Coinc.SphRecons.DistSource(k,:);
                minDistSource(NbCoincFinal)=StructT.Coinc.SphRecons.minDistSource(k);
                Chi2DelaySp(NbCoincFinal)=StructT.Coinc.SphRecons.Chi2Delay(k);
                SlopeDelaySp(NbCoincFinal)=StructT.Coinc.SphRecons.SlopeDelay(k);
                
                if DelayCorrFlag
                    % Coinc.DelayCorr.PlanRecons
                    FlagCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Flag(k);
                    ThetaCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Theta(k);
                    dThetaCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.dTheta(k);
                    PhiCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Phi(k);
                    dPhiCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.dPhi(k);
                    Chi2Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Chi2(k);
                    SignifCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Signif(k);
                    Chi2DelayCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.Chi2Delay(k);
                    SlopeDelayCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.PlanRecons.SlopeDelay(k);
                
                    % Coinc.DelayCorr.SphRecons
                    FlagSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Flag(k);
                    RhoCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Rho(k);
                    ThetaSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Theta(k);
                    PhiSpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Phi(k);
                    X0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.X0(k);
                    Y0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Y0(k);
                    Z0Corr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Z0(k);
                    DistSourceCorr(NbCoincFinal,indAnt)=StructT.Coinc.DelayCorrRecons.SphRecons.DistSource(k,:);
                    Chi2DelaySpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.Chi2Delay(k);
                    SlopeDelaySpCorr(NbCoincFinal)=StructT.Coinc.DelayCorrRecons.SphRecons.SlopeDelay(k);
                end
                
                % Coinc.ShowerRecons
                for m=1:length(Ant)
                    ShowerSignals(NbCoincFinal,indAnt(m))=StructT.Coinc.ShowerRecons.ShowerSignals(k,m);
                end;
                if isfield(StructT.Coinc.ShowerRecons,'XCore')
                    XCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.XCore(k);
                    YCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.YCore(k);
                    ZCore(NbCoincFinal)=StructT.Coinc.ShowerRecons.ZCore(k);
                    AxisAmp(NbCoincFinal)=StructT.Coinc.ShowerRecons.AxisAmp(k);
                    Lambda(NbCoincFinal)=StructT.Coinc.ShowerRecons.Lambda(k);
                end;
            end;
        end;
        CoincMetaIter=CoincMetaIter+CoincMeta;        
    end;
       
    TotalCoinc=TotalCoinc+CoincMetaIter;
    for k=1:length(Ant)
        EvtAntTotal(indAnt(k))=EvtAntTotal(indAnt(k))+StructT.Setup.Det(k).Evt;
    end;
    
end;
            
% Save dst file

Struct.Setup=Setup;
Struct.Setup.TotalCoinc=NbCoincFinal;
Struct.Coinc.Mult=Mult(1:NbCoincFinal)';
Struct.Coinc.MultAnt=MultAnt(1:NbCoincFinal)';
Struct.Coinc.IdCoinc=IdCoinc(1:NbCoincFinal)';
Struct.Coinc.IsShower=IsShower(1:NbCoincFinal)';

Struct.Coinc.Det.Status=Status(1:NbCoincFinal,:);
Struct.Coinc.Det.Evt=Evt(1:NbCoincFinal,:);
Struct.Coinc.Det.UnixTime=UnixTime(1:NbCoincFinal,:);
Struct.Coinc.Det.TriggerRate=TriggerRate(1:NbCoincFinal,:);
Struct.Coinc.Det.Sigma=Sigma(1:NbCoincFinal,:);
Struct.Coinc.Det.MinRaw=MinRaw(1:NbCoincFinal,:);
Struct.Coinc.Det.MaxRaw=MaxRaw(1:NbCoincFinal,:);
Struct.Coinc.Det.AmpMax=AmpMax(1:NbCoincFinal,:);
Struct.Coinc.Det.TrigTime=TrigTime(1:NbCoincFinal,:);
Struct.Coinc.Det.Gain=Gain(1:NbCoincFinal,:);
Struct.Coinc.Det.TrigCor=TrigCor(1:NbCoincFinal,:);
Struct.Coinc.Det.CoefCor=CoefCor(1:NbCoincFinal,:);
Struct.Coinc.Det.GainPSD=GainPSD(1:NbCoincFinal);
Struct.Coinc.Det.CalibratedAmp1=CalibratedAmp1(1:NbCoincFinal,:);
Struct.Coinc.Det.CalibratedAmp2=CalibratedAmp2(1:NbCoincFinal,:);
Struct.Coinc.Det.Data3Scints=Data3Scints(1:NbCoincFinal,:);
                
Struct.Coinc.PlanRecons.Radio.Flag=FlagR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.Theta=ThetaR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.dTheta=dThetaR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.Phi=PhiR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.dPhi=dPhiR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.Chi2=Chi2R(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.Signif=SignifR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.Chi2Delay=Chi2DelayR(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Radio.SlopeDelay=SlopeDelayR(1:NbCoincFinal);
                
Struct.Coinc.PlanRecons.Hybrid.Flag=FlagH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.Theta=ThetaH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.dTheta=dThetaH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.Phi=PhiH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.dPhi=dPhiH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.Chi2=Chi2H(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.Signif=SignifH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.Chi2Delay=Chi2DelayH(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Hybrid.SlopeDelay=SlopeDelayH(1:NbCoincFinal);
                
Struct.Coinc.PlanRecons.Sci.Flag=FlagS(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Sci.Theta=ThetaS(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Sci.dTheta=dThetaS(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Sci.Phi=PhiS(1:NbCoincFinal);
Struct.Coinc.PlanRecons.Sci.dPhi=dPhiS(1:NbCoincFinal);
                
Struct.Coinc.SphRecons.Flag=FlagSp(1:NbCoincFinal);
Struct.Coinc.SphRecons.Rho=Rho(1:NbCoincFinal);
Struct.Coinc.SphRecons.Theta=ThetaSp(1:NbCoincFinal);
Struct.Coinc.SphRecons.Phi=PhiSp(1:NbCoincFinal);
Struct.Coinc.SphRecons.X0=X0(1:NbCoincFinal);
Struct.Coinc.SphRecons.Y0=Y0(1:NbCoincFinal);
Struct.Coinc.SphRecons.Z0=Z0(1:NbCoincFinal);
Struct.Coinc.SphRecons.DistSource=DistSource(1:NbCoincFinal,:);
%Struct.Coinc.SphRecons.minDistSource=minDistSource(1:NbCoincFinal);
Struct.Coinc.SphRecons.Chi2Delay=Chi2DelaySp(1:NbCoincFinal);
Struct.Coinc.SphRecons.SlopeDelay=SlopeDelaySp(1:NbCoincFinal);

if DelayCorrFlag
    Struct.Coinc.DelayCorrRecons.PlanRecons.Flag=FlagCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.Theta=ThetaCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.dTheta=dThetaCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.Phi=PhiCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.dPhi=dPhiCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.Chi2=Chi2Corr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.Signif=SignifCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.Chi2Delay=Chi2DelayCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.PlanRecons.SlopeDelay=SlopeDelayCorr(1:NbCoincFinal);
                
    Struct.Coinc.DelayCorrRecons.SphRecons.Flag=FlagSpCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.Rho=RhoCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.Theta=ThetaSpCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.Phi=PhiSpCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.X0=X0Corr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.Y0=Y0Corr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.Z0=Z0Corr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.DistSource=DistSourceCorr(1:NbCoincFinal,:);
    Struct.Coinc.DelayCorrRecons.SphRecons.Chi2Delay=Chi2DelaySpCorr(1:NbCoincFinal);
    Struct.Coinc.DelayCorrRecons.SphRecons.SlopeDelay=SlopeDelaySpCorr(1:NbCoincFinal);
end;

Struct.Coinc.ShowerRecons.ShowerSignals=ShowerSignals(1:NbCoincFinal,:);
if isfield(StructT.Coinc.ShowerRecons,'XCore')
    Struct.Coinc.ShowerRecons.XCore=XCore(1:NbCoincFinal);
    Struct.Coinc.ShowerRecons.YCore=YCore(1:NbCoincFinal);
    Struct.Coinc.ShowerRecons.ZCore=ZCore(1:NbCoincFinal);
    Struct.Coinc.ShowerRecons.AxisAmp=AxisAmp(1:NbCoincFinal);
    Struct.Coinc.ShowerRecons.Lambda=Lambda(1:NbCoincFinal);
end;

Struct.Coinc.InfosSignal.NbBox=NbBoxSave(1:NbCoincFinal,:);
Struct.Coinc.InfosSignal.TotalToT=TotalToTSave(1:NbCoincFinal,:);
Struct.Coinc.InfosSignal.BoxToT=BoxToTSave(1:NbCoincFinal,:);
Struct.Coinc.InfosSignal.TimeDiff=TimeDiffSave(1:NbCoincFinal);
                
filename=[DST_PATH sprintf(dst_filename,nrun,IterFinal)]
save(filename,'Struct')
clear Struct


        
        


