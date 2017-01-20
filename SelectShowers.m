function [ output_args ] = SelectShowers( nrun )
% Selection of showers and displays
% OMH 12/0

SharedGlobals;
mtyp= ['+bla';'.blu';'dmag';'pred';'hgre';'oyel'];
htyp = ['sb';'sk';'sr'];
linf = 4;
lsup = 14;
step = 2;
    
%% Load dst
dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
dst = load(dstname);
dst.Struct = Dist2Source(dst.Struct);

ncoincs = dst.Struct.Setup.TotalCoinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
DetectorType = [dst.Struct.Setup.Det.isScint];
nDets = length(Detectors);
nScints = sum(DetectorType);
nAnts = nDets-nScints;
Events = [DetStruct.Evt];

RunSetup = dst.Struct.Setup;
DetPosX = [RunSetup.Det.X];
DetPosY = [RunSetup.Det.Y];
DetPosZ = [RunSetup.Det.Z];

% Recons
CoincStruct = dst.Struct.Coinc;
Id =  CoincStruct.IdCoinc;
PlanStruct = CoincStruct.PlanRecons.Radio;
SphStruct = CoincStruct.SphRecons;
mult = PlanStruct.L;
time = SphStruct.T/FSAMPLING/60;  % minutes
idp = CoincStruct.IdCoinc;
thetap = PlanStruct.Theta;
phip = PlanStruct.Phi;
chi2p = PlanStruct.Chi2/(mult-2);
thetas = SphStruct.Theta;
phis = SphStruct.Phi;
rhos = SphStruct.Rho;
xs = SphStruct.X0;
ys = SphStruct.Y0;
zs = SphStruct.Z0;
r = SphStruct.minDistSource;

% ShowerRecons
ShowerStruct = CoincStruct.ShowerRecons;
isShower = CoincStruct.IsShower;
iShower = find(isShower == 1);
nShowers = length(iShower);
lambda = ShowerStruct.Lambda;

minDistCore = zeros(1,ncoincs);
for i = 1:ncoincs
    if isShower(i)==1
        distCore = ([ShowerStruct.XCore(i)-DetPosX ShowerStruct.YCore(i)-DetPosY ShowerStruct.ZCore(i)-DetPosZ] );
        minDistCore(i) = min(norm(distCore));
    end
end
size(r)
size(mult)
sel = mult>6 & minDistCore < 200 & thetap < 65 & thetap-thetas<5 & phip-phis<5 & r'>500; 

for i = 1:length(sel)
    AnaCoinc(nrun,Id(sel(i)),0)
    pause
    close all
end



end

