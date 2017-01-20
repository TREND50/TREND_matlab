function [Struct]=EmptyStructure(filename,Struct,nrun,NbIter)

SharedGlobals;

% Initialisation empty structures for PreProcess and Dst file
StructT=Struct;
clear Struct;

if strcmp(filename,sprintf(coinc_filename,nrun,NbIter))
    
    Struct.Setup=StructT.Setup;
    Struct.Setup.TotalCoinc=0;
    Struct.Setup.InfosRun.GlobalCoincRateQuickReject=0;
    Struct.Setup.InfosRun.DetCoincRateQuickReject=0;
    Struct.Setup.InfosRun.GlobalCoincRateRaw=0;
    Struct.Setup.InfosRun.DetCoincRateRaw=0;
    
    % Structure Coinc
    Coinc.Mult=0;
    Coinc.MultAnt=0;
    Coinc.MultSci=0;
    Coinc.IdCoinc=0;
    
    Coinc.Det.Tag=0;
    Coinc.Det.Status=0;
    Coinc.Det.Id=0;
    Coinc.Det.Evt=0;
    Coinc.Det.Time=0;
    Coinc.Det.UnixTime=0;
    Coinc.Det.TriggerRate=0;
    
    Struct.Coinc=Coinc;
    
    filename = [DST_PATH sprintf(coinc_filename,nrun,NbIter)];
    save(filename,'Struct');
    
elseif strcmp(filename,sprintf(preprocess_filename,nrun,NbIter))

    Struct.Setup=StructT.Setup;
    Struct.Setup.TotalCoinc=0;
    Struct.Setup.InfosRun.GlobalCoincRateQuickReject=0;
    Struct.Setup.InfosRun.DetCoincRateQuickReject=0;
    Struct.Setup.InfosRun.GlobalCoincRateRaw=0;
    Struct.Setup.InfosRun.DetCoincRateRaw=0;
    
    % Structure Coinc
    Coinc.Mult=0;
    Coinc.MultAnt=0;
    Coinc.MultSci=0;
    Coinc.IdCoinc=0;

    Coinc.Det.Tag=0;
    Coinc.Det.Status=0;
    Coinc.Det.Id=0;
    Coinc.Det.Evt=0;
    Coinc.Det.Time=0;
    Coinc.Det.UnixTime=0;
    Coinc.Det.TriggerRate=0;
    Coinc.Det.Sigma=0;
    Coinc.Det.Mu=0;
    Coinc.Det.MinRaw=0;
    Coinc.Det.MaxRaw=0;
    Coinc.Det.Sat=0;
    Coinc.Det.TimeSave=0;

    Coinc.Reject.ConsCoinc=0;
    Coinc.Reject.RawFilter=0;
    Coinc.Reject.NewMult=0;
    Coinc.Reject.CaracEvt=0;
    Coinc.Reject.NoBox=0;
    
    Struct.Coinc=Coinc;
    
    filename = [DST_PATH sprintf(preprocess_filename,nrun,NbIter)];
    save(filename,'Struct');
    
elseif strcmp(filename,sprintf(dst_filename,nrun,NbIter))
    
    Struct.Setup=StructT.Setup;
    Struct.Setup.TotalCoinc=0;
    
    % Structure Coinc
    Coinc.Mult=0;
    Coinc.MultAnt=0;
    Coinc.MultSci=0;
    Coinc.IdCoinc=0;
    Coinc.IsShower=0;
    
    Coinc.Det.Tag=0;
    Coinc.Det.Status=0;
    Coinc.Det.Id=0;
    Coinc.Det.Evt=0;
    Coinc.Det.Time=0;
    Coinc.Det.UnixTime=0;
    Coinc.Det.TriggerRate=0;
    Coinc.Det.Sigma=0;
    Coinc.Det.Mu=0;
    Coinc.Det.MinRaw=0;
    Coinc.Det.MaxRaw=0;
    Coinc.Det.Sat=0;
    Coinc.Det.AmpMax=0;
    Coinc.Det.TrigTime=0;
    Coinc.Det.Gain=0;
    Coinc.Det.TrigCor=0;
    Coinc.Det.CoefCor=0;
    Coinc.Det.GainPSD=0;
    Coinc.Det.CalibratedAmp1=0;
    Coinc.Det.CalibratedAmp2=0;
    Coinc.Det.Data3Scints=0;
    
    Coinc.PlanRecons.Radio.Flag=0;
    Coinc.PlanRecons.Radio.Theta=0;
    Coinc.PlanRecons.Radio.dTheta=0;
    Coinc.PlanRecons.Radio.Phi=0;
    Coinc.PlanRecons.Radio.dPhi=0;
    Coinc.PlanRecons.Radio.Chi2=0;
    Coinc.PlanRecons.Radio.Signif=0;
    Coinc.PlanRecons.Radio.Chi2Delay=0;
    Coinc.PlanRecons.Radio.SlopeDelay=0;
    
    Coinc.PlanRecons.Hybrid.Flag=0;
    Coinc.PlanRecons.Hybrid.Theta=0;
    Coinc.PlanRecons.Hybrid.dTheta=0;
    Coinc.PlanRecons.Hybrid.Phi=0;
    Coinc.PlanRecons.Hybrid.dPhi=0;
    Coinc.PlanRecons.Hybrid.Chi2=0;
    Coinc.PlanRecons.Hybrid.Signif=0;
    Coinc.PlanRecons.Hybrid.Chi2Delay=0;
    Coinc.PlanRecons.Hybrid.SlopeDelay=0;
    
    Coinc.PlanRecons.Sci.Flag=0;
    Coinc.PlanRecons.Sci.Theta=0;
    Coinc.PlanRecons.Sci.dTheta=0;
    Coinc.PlanRecons.Sci.Phi=0;
    Coinc.PlanRecons.Sci.dPhi=0;
    
    Coinc.SphRecons.Flag=0;
    Coinc.SphRecons.Rho=0;
    Coinc.SphRecons.Theta=0;
    Coinc.SphRecons.Phi=0;
    Coinc.SphRecons.X0=0;
    Coinc.SphRecons.Y0=0;
    Coinc.SphRecons.Z0=0;
    Coinc.SphRecons.DistSource=0;
    Coinc.SphRecons.minDistSource=0;
    Coinc.SphRecons.Chi2Delay=0;
    Coinc.SphRecons.SlopeDelay=0;
    
    Coinc.ShowerRecons.ShowerSignals=0;
    
    Coinc.ShowerRecons.RawRecons.XCore=0;
    Coinc.ShowerRecons.RawRecons.YCore=0;
    Coinc.ShowerRecons.RawRecons.ZCore=0;
    Coinc.ShowerRecons.RawRecons.AxisAmp=0;
    Coinc.ShowerRecons.RawRecons.Lambda=0; 
    
    Coinc.ShowerRecons.CalRecons.XCore=0;
    Coinc.ShowerRecons.CalRecons.YCore=0;
    Coinc.ShowerRecons.CalRecons.ZCore=0;
    Coinc.ShowerRecons.CalRecons.AxisAmp=0;
    Coinc.ShowerRecons.CalRecons.Lambda=0; 
    
    
    Coinc.DelayCorrRecons.IsShower=0;
    
    Coinc.DelayCorrRecons.PlanRecons.Radio.Flag=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.Theta=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.dTheta=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.Phi=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.dPhi=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.Chi2=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.Signif=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.Chi2Delay=0;
    Coinc.DelayCorrRecons.PlanRecons.Radio.SlopeDelay=0;
    
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Flag=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Theta=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.dTheta=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Phi=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.dPhi=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Chi2=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Signif=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.Chi2Delay=0;
    Coinc.DelayCorrRecons.PlanRecons.Hybrid.SlopeDelay=0;
    
    Coinc.DelayCorrRecons.SphRecons.Flag=0;
    Coinc.DelayCorrRecons.SphRecons.Rho=0;
    Coinc.DelayCorrRecons.SphRecons.Theta=0;
    Coinc.DelayCorrRecons.SphRecons.Phi=0;
    Coinc.DelayCorrRecons.SphRecons.X0=0;
    Coinc.DelayCorrRecons.SphRecons.Y0=0;
    Coinc.DelayCorrRecons.SphRecons.Z0=0;
    Coinc.DelayCorrRecons.SphRecons.DistSource=0;
    Coinc.DelayCorrRecons.SphRecons.minDistSource=0;
    Coinc.DelayCorrRecons.SphRecons.Chi2Delay=0;
    Coinc.DelayCorrRecons.SphRecons.SlopeDelay=0;
    
    Coinc.DelayCorrRecons.ShowerRecons.ShowerSignals=0;
    
    Coinc.DelayCorrRecons.ShowerRecons.RawRecons.XCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.RawRecons.YCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.RawRecons.ZCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.RawRecons.AxisAmp=0;
    Coinc.DelayCorrRecons.ShowerRecons.RawRecons.Lambda=0; 
    
    Coinc.DelayCorrRecons.ShowerRecons.CalRecons.XCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.CalRecons.YCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.CalRecons.ZCore=0;
    Coinc.DelayCorrRecons.ShowerRecons.CalRecons.AxisAmp=0;
    Coinc.DelayCorrRecons.ShowerRecons.CalRecons.Lambda=0; 
    
    Struct.Coinc=Coinc;
    
    filename = [DST_PATH sprintf(dst_filename,nrun,NbIter)];
    save(filename,'Struct');
    
end;
    
