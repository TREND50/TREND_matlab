function [Struct]=CleanFilt(Struct,slice)
% Cut on nboxes and pulse duration
% OMH 06/12/2011

SharedGlobals;

nrun = Struct.Setup.Run;
PulseStruct = Struct.Coinc.InfosSignal;
nbboxes = PulseStruct.NbBox;
ttot = PulseStruct.TotalToT;
tcen = PulseStruct.BoxToT;

hist(ttot(ttot>0))
NbCoinc = Struct.Setup.TotalCoinc;

lsel = [];

for i=1:NbCoinc
    in = [Struct.Coinc.Det.Status(i,:) >-1];
    if isempty(find(nbboxes(i,in)==0)) & isempty(find(nbboxes(i,in)>3)) & isempty(find(ttot(i,in)>0.6e-6)) & isempty(find(tcen(i,in)>0.6e-6))        
        lsel(end+1) = i;  % Tag all coincs with at list 
    end
end

Struct.Setup.TotalCoinc = length(lsel);

% Coinc struc
Struct.Coinc.Mult = Struct.Coinc.Mult(lsel);
Struct.Coinc.MultAnt = Struct.Coinc.MultAnt(lsel);
Struct.Coinc.IdCoinc = Struct.Coinc.IdCoinc(lsel);
Struct.Coinc.IsShower = Struct.Coinc.IsShower(lsel);
% Det struct
Struct.Coinc.Det.Status = Struct.Coinc.Det.Status(lsel,:);
Struct.Coinc.Det.Evt = Struct.Coinc.Det.Evt(lsel,:);
Struct.Coinc.Det.UnixTime = Struct.Coinc.Det.UnixTime(lsel,:);
Struct.Coinc.Det.TriggerRate = Struct.Coinc.Det.TriggerRate(lsel,:);
Struct.Coinc.Det.Sigma = Struct.Coinc.Det.Sigma(lsel,:);
Struct.Coinc.Det.MinRaw = Struct.Coinc.Det.MinRaw(lsel,:);
Struct.Coinc.Det.MaxRaw = Struct.Coinc.Det.MaxRaw(lsel,:);
Struct.Coinc.Det.AmpMax = Struct.Coinc.Det.AmpMax(lsel,:);
Struct.Coinc.Det.TrigTime = Struct.Coinc.Det.TrigTime(lsel,:);
Struct.Coinc.Det.Gain = Struct.Coinc.Det.Gain(lsel,:);
Struct.Coinc.Det.TrigCor = Struct.Coinc.Det.TrigCor(lsel,:);
Struct.Coinc.Det.CoefCor = Struct.Coinc.Det.CoefCor(lsel,:);
Struct.Coinc.Det.CalibratedAmp1 = Struct.Coinc.Det.CalibratedAmp1(lsel,:);
Struct.Coinc.Det.CalibratedAmp2 = Struct.Coinc.Det.CalibratedAmp2(lsel,:);

% PlanRecons struct
Struct.Coinc.PlanRecons.Radio.Flag = Struct.Coinc.PlanRecons.Radio.Flag(lsel);
Struct.Coinc.PlanRecons.Radio.Theta = Struct.Coinc.PlanRecons.Radio.Theta(lsel);
Struct.Coinc.PlanRecons.Radio.dTheta = Struct.Coinc.PlanRecons.Radio.dTheta(lsel);
Struct.Coinc.PlanRecons.Radio.Phi = Struct.Coinc.PlanRecons.Radio.Phi(lsel);
Struct.Coinc.PlanRecons.Radio.dPhi = Struct.Coinc.PlanRecons.Radio.dPhi(lsel);
Struct.Coinc.PlanRecons.Radio.Chi2 = Struct.Coinc.PlanRecons.Radio.Chi2(lsel);
Struct.Coinc.PlanRecons.Radio.Signif = Struct.Coinc.PlanRecons.Radio.Signif(lsel);
Struct.Coinc.PlanRecons.Radio.Chi2Delay = Struct.Coinc.PlanRecons.Radio.Chi2Delay(lsel);
Struct.Coinc.PlanRecons.Radio.SlopeDelay = Struct.Coinc.PlanRecons.Radio.SlopeDelay(lsel);
%
Struct.Coinc.Planecons.Hybrid.Flag = Struct.Coinc.PlanRecons.Hybrid.Flag(lsel);
Struct.Coinc.PlanRecons.Hybrid.Theta = Struct.Coinc.PlanRecons.Hybrid.Theta(lsel);
Struct.Coinc.PlanRecons.Hybrid.dTheta = Struct.Coinc.PlanRecons.Hybrid.dTheta(lsel);
Struct.Coinc.PlanRecons.Hybrid.Phi = Struct.Coinc.PlanRecons.Hybrid.Phi(lsel);
Struct.Coinc.PlanRecons.Hybrid.dPhi = Struct.Coinc.PlanRecons.Hybrid.dPhi(lsel);
Struct.Coinc.PlanRecons.Hybrid.Chi2 = Struct.Coinc.PlanRecons.Hybrid.Chi2(lsel);
Struct.Coinc.PlanRecons.Hybrid.Signif = Struct.Coinc.PlanRecons.Hybrid.Signif(lsel);
Struct.Coinc.PlanRecons.Hybrid.Chi2Delay = Struct.Coinc.PlanRecons.Hybrid.Chi2Delay(lsel);
Struct.Coinc.PlanRecons.Hybrid.SlopeDelay = Struct.Coinc.PlanRecons.Hybrid.SlopeDelay(lsel);
%
Struct.Coinc.Planecons.Sci.Flag = Struct.Coinc.PlanRecons.Sci.Flag(lsel);
Struct.Coinc.PlanRecons.Sci.Theta = Struct.Coinc.PlanRecons.Sci.Theta(lsel);
Struct.Coinc.PlanRecons.Sci.dTheta = Struct.Coinc.PlanRecons.Sci.dTheta(lsel);
Struct.Coinc.PlanRecons.Sci.Phi = Struct.Coinc.PlanRecons.Sci.Phi(lsel);
Struct.Coinc.PlanRecons.Sci.dPhi = Struct.Coinc.PlanRecons.Sci.dPhi(lsel);
% SphRecons struct
Struct.Coinc.SphRecons.Flag = Struct.Coinc.SphRecons.Flag(lsel);
Struct.Coinc.SphRecons.Rho = Struct.Coinc.SphRecons.Rho(lsel);
Struct.Coinc.SphRecons.Theta = Struct.Coinc.SphRecons.Theta(lsel);
Struct.Coinc.SphRecons.Phi = Struct.Coinc.SphRecons.Phi(lsel);
Struct.Coinc.SphRecons.X0 = Struct.Coinc.SphRecons.X0(lsel);
Struct.Coinc.SphRecons.Y0 = Struct.Coinc.SphRecons.Y0(lsel);
Struct.Coinc.SphRecons.Z0 = Struct.Coinc.SphRecons.Z0(lsel);
Struct.Coinc.SphRecons.DistSource = Struct.Coinc.SphRecons.DistSource(lsel,:);
Struct.Coinc.SphRecons.Chi2Delay = Struct.Coinc.SphRecons.Chi2Delay(lsel);
Struct.Coinc.SphRecons.SlopeDelay = Struct.Coinc.SphRecons.SlopeDelay(lsel);
% ShowerRecons struct
Struct.Coinc.ShowerRecons.XCore = Struct.Coinc.ShowerRecons.XCore(lsel);
Struct.Coinc.ShowerRecons.YCore = Struct.Coinc.ShowerRecons.YCore(lsel);
Struct.Coinc.ShowerRecons.ZCore = Struct.Coinc.ShowerRecons.ZCore(lsel);
Struct.Coinc.ShowerRecons.AxisAmp = Struct.Coinc.ShowerRecons.AxisAmp(lsel);
Struct.Coinc.ShowerRecons.Lambda = Struct.Coinc.ShowerRecons.Lambda(lsel);
% InfosSignal struct
Struct.Coinc.InfosSignal.NbBox = Struct.Coinc.InfosSignal.NbBox(lsel,:);
Struct.Coinc.InfosSignal.TotalToT = Struct.Coinc.InfosSignal.TotalToT(lsel,:);
Struct.Coinc.InfosSignal.BoxToT = Struct.Coinc.InfosSignal.BoxToT(lsel,:);
Struct.Coinc.InfosSignal.TimeDiff = Struct.Coinc.InfosSignal.TimeDiff(lsel);
% Cells
for i=1:length(lsel)
    NewGainPSD{i} = Struct.Coinc.Det.GainPSD{lsel(i)};    
    for j = 1:size(Struct.Coinc.ShowerRecons.ShowerSignals,2)
        NewShowerSignals{i,j} = Struct.Coinc.ShowerRecons.ShowerSignals{lsel(i),j};
        NewData3Scints{i,j} = Struct.Coinc.Det.Data3Scints{lsel(i),j};    
    end
end
Struct.Coinc.Det.GainPSD = NewGainPSD;
Struct.Coinc.ShowerRecons.ShowerSignals = NewShowerSignals;
Struct.Coinc.Det.Data3Scints = NewData3Scints;

filename =['D:/dst/cleanfilt/' sprintf(dst_filename,nrun,slice)];
save(filename,'Struct');
display(sprintf('DST %s now saved to file.',filename));
display(sprintf('%d coincs saved out of %d initially.',length(lsel),NbCoinc));