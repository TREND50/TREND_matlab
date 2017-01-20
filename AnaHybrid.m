function [ ] = AnaHybrid( nrun )
%Analyse hybrid dst
% OMH 08/08/2011

SharedGlobals;

%% Load dst
dstname = [DST_PATH 'hybrid/' sprintf(dsthyb_filename,nrun,1)];
if ~exist(dstname)
    disp(sprintf('dst %s not found',dstname))
return
end
dst = load(dstname);

Struct = dst.Struct;
Struct = Dist2Source(Struct);

%% Get DST variables
Detectors = [Struct.Setup.Det.Name];
DetectorType = [Struct.Setup.Det.isScint];
SetupStruct = Struct.Setup;
CoincStruct = Struct.Coinc;
Evt = CoincStruct.Det.Evt;
ReconsHybStruct = CoincStruct.PlanRecons.Hybrid;
ReconsSciStruct = CoincStruct.PlanRecons.Sci;
ReconsAntStruct = CoincStruct.PlanRecons.Radio;
tag = [CoincStruct.Det.Tag];
multAnt = CoincStruct.MultAnt;
multSci = CoincStruct.MultSci;

isHybrid = find(ReconsHybStruct.Flag==1);

theta_ant = ReconsAntStruct.Theta(isHybrid);
phi_ant = ReconsAntStruct.Phi(isHybrid);
theta_hyb = ReconsHybStruct.Theta(isHybrid);
phi_hyb = ReconsHybStruct.Phi(isHybrid);
theta_sci = ReconsSciStruct.Theta(isHybrid);
phi_sci = ReconsSciStruct.Phi(isHybrid);
tstart = SetupStruct.InfosRun.TimeStart(1);
sec = max(CoincStruct.Det.UnixTime(isHybrid,:),[],2);
thyb = (sec-tstart)/TriggerTimeSpan;
nhybs = length(isHybrid);
icoinc_hyb = CoincStruct.IdCoinc(isHybrid);

disp(sprintf('R%d : %d hybrids',nrun,nhybs));
for i = 1:nhybs
    iant_hyb = find(tag(isHybrid(i),:)==1 & DetectorType == 0);
    isci_hyb = find(tag(isHybrid(i),:)==1 & DetectorType == 1);
    na = length(iant_hyb);
    ns = length(isci_hyb);
    disp(sprintf('Coinc %d : %d antennas + %d scints',icoinc_hyb(i),na,ns))
    Detectors(iant_hyb)
    Evt(isHybrid(i),iant_hyb)
    Detectors(isci_hyb)
    Evt(isHybrid(i),isci_hyb)

    %% Trigger rates at hybrid time
    disp(sprintf('Unix second = %d, time since run start = %3.1f mins',sec(i),thyb(i))); %Assuming TriggerTimeSpan = 60;  
    % Get antenna data
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    dstr = load(dstname);
    coincrate = dstr.Struct.Setup.InfosRun.GlobalCoincRateRaw;
    trigrate = dstr.Struct.Setup.InfosRun.TrigRate;
    tagr = [dstr.Struct.Coinc.Det.Tag];
    Detectorr = [dstr.Struct.Setup.Det.Name];
    i139 = find(tagr(:,Detectorr==139)==1);
    i175 = find(tagr(:,Detectorr==175)==1);
    i176 = find(tagr(:,Detectorr==176)==1);
    i1=intersect(i139,i175);
    i2=intersect(i139,i176);
    i3= intersect(i176,i176);
    triples = length(intersect(i139,i2));
    doubles = length(i1)+length(i2)+length(i3)-triples;
    total_scints = doubles+triples;
    disp 'Det      Trig rate '
    for j = 1:multAnt(isHybrid(i))
        disp(sprintf('%d     %6.3f Hz ',Detectors(iant_hyb(j)),trigrate(round(thyb(i)),iant_hyb(j))))
    end
    disp(sprintf('Radio coinc rate at %d mins = %6.3f Hz',round(thyb(i)),coincrate(round(thyb(i)))))        
    for j = 1:multSci(isHybrid(i))
        disp(sprintf('%d     %6.3f Hz ',Detectors(isci_hyb(j)),trigrate(round(thyb(i)),isci_hyb(j))))
    end
    disp(sprintf('Scint coinc rate over run = %d/%3.1f = %1.3f /min',total_scints,length(coincrate),total_scints/length(coincrate)))        
    
    figure(12)
    subplot(2,1,1)
    plot(coincrate)
    xlim([0 length(coincrate)])
    xlabel('Time [min]')
    ylabel('Trigger Rate [Hz]')
    subplot(2,1,2)
    for j = 1:multAnt(isHybrid(i))
        plot(trigrate(:,iant_hyb(j)),'k','LineWidth',2)
        xlim([0 length(coincrate)])
        xlabel('Time [min]')
        ylabel('Trigger Rate [Hz]')
        line([thyb(i) thyb(i)],[0 max(trigrate(:,iant_hyb(j)))])
        evt = [dstr.Struct.Setup.Det.Evt];
        %disp(sprintf('Detector %d , total event number = %d (%3.1f)', Detectors(iant_hyb(j)),evt(iant_hyb(j)),mean(trigrate(:,iant_hyb(j)))*60*length(coincrate)))
        hold on
    end
    clear dstr
    AnaRecons(nrun,round(thyb(i)))
    
    %% Reconstruction
    if na>3
        disp(sprintf('Antennes: Theta = %3.1f, Phi = %3.1f',theta_ant(i),phi_ant(i)))
    end
    if ns==3
        disp(sprintf('Scintillators: Theta = %3.1f, Phi = %3.1f',theta_sci(i),phi_sci(i)))
    end
    disp(sprintf('Hybrid: Theta = %3.1f, Phi = %3.1f',theta_hyb(i),phi_hyb(i)))
    if na>3
        AnaCoinc(nrun,icoinc_hyb(i),1)
    end
    pause
    close all
end