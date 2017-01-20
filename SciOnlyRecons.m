function [Struct]=SciOnlyRecons(Struct, AnalysisType)
% Arrival direction reconstruction for scintillators only
% TS

SharedGlobals;

indEvt=find([Struct.Coinc.MultSci]==3);
ncoinc=Struct.Setup.TotalCoinc;
nrun = Struct.Setup.Run;

display(sprintf('%d coincs with 3 scints...',length(indEvt)))

indsci=find([Struct.Setup.Det.isScint]==1);
posX=[Struct.Setup.Det(indsci).Y];
posY=-1.*[Struct.Setup.Det(indsci).X];
posZ=[Struct.Setup.Det(indsci).Z];
time=[Struct.Coinc.Det.TrigTime];
Evt = [Struct.Coinc.Det.Evt];
Det = [Struct.Setup.Det.Name];
flag=zeros(1,ncoinc);
theta=-1*ones(1,ncoinc);
phi=-1*ones(1,ncoinc);
dtheta=-1*ones(1,ncoinc);
dphi=-1*ones(1,ncoinc);
Data3Scints = cell(ncoinc,length(Det));

for i=1:length(indEvt)
    
    if AnalysisType ==2
        % Get scints raw data and write it to structure
        for j=1:length(indsci)
            fdc = OpenFileData( nrun, Det(indsci(j)));
            if fdc>0
                ev = Evt(indEvt(i),indsci(j));
                fseek(fdc,(ev-1)*ibuffs,'bof'); % Skip events
                Data3Scints{indEvt(i),indsci(j)} = fread(fdc,ibuffs,'*uint8');
            end
        end
    end
    
    param=zeros(3,5);
    param(:,1) = indsci;
    param(:,2) = posX;  %SN
    param(:,3) = posY; %EW = -WE
    param(:,4) = posZ;
    param(:,5) = time(indEvt(i),indsci);

    [thetat phit] = fitplan_3(param);
    [dthetat dphit] = ComputeScintDirError(Struct,thetat,phit);
    if isreal(phit)
        flag(indEvt(i))=1;
        if thetat>90
            thetat = 180-thetat;
        end
        theta(indEvt(i))=thetat;
        dtheta(indEvt(i))=dthetat;
        
        phit = mod(phit,360);
        phi(indEvt(i))=phit;
        dphi(indEvt(i))=dphit;
    end;
    
end;
Struct.Coinc.PlanRecons.Sci.Flag=flag;
Struct.Coinc.PlanRecons.Sci.Theta=theta;
Struct.Coinc.PlanRecons.Sci.dTheta=dtheta;
Struct.Coinc.PlanRecons.Sci.Phi=phi;   
Struct.Coinc.PlanRecons.Sci.dPhi=dphi;   
Struct.Coinc.Det.Data3Scints = Data3Scints;    

