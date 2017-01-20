function [] = BuildLightDST(nrun)
% Remove signals from DSTs to make them lighter
% OMH 11/06/2013

SharedGlobals;


%% Get number of ub dsts
stopflag=0;
nbiter=1;
while stopflag==0
    filename = [DST_PATH sprintf('dst%d_%d.mat',nrun,nbiter)];
    fd=fopen(filename);
    if fd~=-1
        nbiter=nbiter+1;
        fclose(fd);
    else
        stopflag=1;
    end;
end;
nbiter=nbiter-1;

if nbiter==0
    display(sprintf('No dst found for run %d',nrun))
end;

%% Loop on sub dsts
for j=1:nbiter

    filename = [DST_PATH sprintf('dst%d_%d.mat',nrun,j)];
    display(sprintf('Loading dst%d_%d.mat ...',nrun,j))
    dst = load(filename);
    disp 'Done.'
    dst.Struct.Coinc.Det.Data3Scints = {};
    dst.Struct.Coinc.Det.CalibratedAmp1 = [];    
    dst.Struct.Coinc.Det.CalibratedAmp2 = [];
    dst.Struct.Coinc.Det.CoefCor = [];
    dst.Struct.Coinc.Det.TriggerRate = [];
    dst.Struct.Coinc.Det.FiltResult = [];
    dst.Struct.Coinc.Det.TrigTime = [];
    %dst.Struct.Coinc.Det.Sigma = [];
    %dst.Struct.Coinc.Det.AmpMax = [];
    dst.Struct.Coinc.Det.MinRaw = [];
    dst.Struct.Coinc.Det.MaxRaw = [];
    
    dst.Struct.Setup.InfosRun.DetCoincRateRaw = [];
    dst.Struct.Setup.InfosRun.TrigRate = [];
    
    
    dst.Struct.Coinc.ShowerRecons = {};
    dst.Struct.Coinc.Candidates = {};
    dst.Struct.Coinc.PlanRecons.Sci = {};
    dst.Struct.Coinc.PlanRecons.Hybrid = {};
    %dst.Struct.Coinc.SphRecons.DistSource = [];
    %dst.Struct.Coinc.DelayCorrRecons = {};
    dst.Struct.Coinc.PlanRecons = {};
    dst.Struct.Coinc.SphRecons = {};
    
    filename = [DST_PATH sprintf('dst%d_%d_light.mat',nrun,j)];
    Struct = dst.Struct;
    save(filename,'Struct');
end