function [Struct]=ReconsLoader(Struct,Type)
% Load C++ reconstruction files and load them into dst
% 22/12/10 OMH

% AnalysisType=0 : radio only (MultAnt>=4, potential sci tag removed)
% AnalysisType=1 : hybrid
% AnalysisType=2 : scintillators

SharedGlobals;

%% Init structures
nrun=Struct.Setup.Run;
ncoinc=Struct.Setup.TotalCoinc;

if Type==0
    filename = [TEXT_PATH sprintf( 'R%d_planerecons.txt', nrun)];
    
    flag=zeros(1,ncoinc);
    timep=-1*ones(1,ncoinc);
    multp=-1*ones(1,ncoinc);
    thetap=-1*ones(1,ncoinc);
    dThetap=-1*ones(1,ncoinc);
    phip=-1*ones(1,ncoinc);
    dPhip=-1*ones(1,ncoinc);
    chi2=-1*ones(1,ncoinc);
    signif=-1*ones(1,ncoinc);
    
    flags=zeros(1,ncoinc);
    times=-1*ones(1,ncoinc);
    mults=-1*ones(1,ncoinc);
    x0=-1*ones(1,ncoinc);
    y0=-1*ones(1,ncoinc);
    z0=-1*ones(1,ncoinc);
    rhos=-1*ones(1,ncoinc);
    thetas=-1*ones(1,ncoinc);
    phis=-1*ones(1,ncoinc);
else
    filename = [TEXT_PATH sprintf( 'R%d_hybrid.txt', nrun)];
    flag=zeros(1,ncoinc);
    timep=-1*ones(1,ncoinc);
    multp=-1*ones(1,ncoinc);
    thetap=-1*ones(1,ncoinc);
    dThetap=-1*ones(1,ncoinc);
    phip=-1*ones(1,ncoinc);
    dPhip=-1*ones(1,ncoinc);
    chi2=-1*ones(1,ncoinc);
    signif=-1*ones(1,ncoinc);
end;

%% Load recons file
disp(sprintf('Opening file %s...',filename)); 
if fopen( filename )>0
    planRes = load( filename );        
else
    disp(sprintf('File %s does not exist.',filename));
    planRecons.Flag = flag;
    planRecons.L = multp;
    planRecons.T = timep;
    planRecons.Theta = thetap;
    planRecons.dTheta = dThetap;
    planRecons.Phi = phip;
    planRecons.dPhi = dPhip;
    planRecons.Chi2 = chi2;
    planRecons.Signif = signif;    
    planRecons.Chi2Delay= chi2; %SL
    planRecons.SlopeDelay= chi2; %SL
    % 
    if Type==0
        Struct.Coinc.PlanRecons.Radio = planRecons;

        sphRecons.Flag = flag;
        sphRecons.L = mults;
        sphRecons.T = times;
        sphRecons.Rho = rhos;
        sphRecons.Theta = thetas;
        sphRecons.Phi = phis;
        sphRecons.X0 = x0;
        sphRecons.Y0 = y0;
        sphRecons.Z0 = z0;

        Struct.Coinc.SphRecons = sphRecons;
    
    elseif Type==1
        Struct.Coinc.PlanRecons.Hybrid = planRecons;
    end;
    return
end  
disp 'Done.'

%%
DetId = [Struct.Coinc.Det.Id];
Tag = [Struct.Coinc.Det.Tag];
TrigTime = [Struct.Coinc.Det.TrigTime];
TrigCor = [Struct.Coinc.Det.TrigCor];
X = [Struct.Setup.Det.X];
Y = [Struct.Setup.Det.Y];
Z = [Struct.Setup.Det.Z];
errort = 1*ErrorTrig/FSAMPLING*1e9; 
chi2ndfp = -1*ones(1,ncoinc);
slopep = -1*ones(1,ncoinc);
chi2ndfs = -1*ones(1,ncoinc);
slopes = -1*ones(1,ncoinc);

for i=1:size(planRes,1)
    
    ind = find(Struct.Coinc.IdCoinc==planRes(i,1));
    flag(ind)=1;
    timep(ind) = planRes(i,2);
    multp(ind) = planRes(i,3);
    thetap(ind) = planRes(i,4);
    dThetap(ind) = planRes(i,5);
    phip(ind) = planRes(i,6);
    dPhip(ind) = planRes(i,7);
    chi2(ind) = planRes(i,8);
    signif(ind) = planRes(i,9);
    %
    thetap(ind) = 180-thetap(ind);
    phip(ind) = phip(ind)+180;
    phip(ind) = mod(phip(ind),360);
    if thetap(ind)>90
        thetap(ind)=180-thetap(ind);
    end;
    %
    tag = Tag(ind,:);
    detId=DetId(ind,tag==1);
    detPos = [X(tag==1)' Y(tag==1)' Z(tag==1)'];
    k= [sind(phip(ind))*sind(thetap(ind)) -cosd(phip(ind))*sind(thetap(ind)) -cosd(thetap(ind))]';  % Warning... With our X-Y convention, cartesian vector coordinates are different from usual 
    plandelays = [ detId' detPos*k ];  % In meters
    plandelays(:,2) = plandelays(:,2) - min(plandelays(:,2)); % Use 1st triggering antenna as ref
    plandelays(:,2) = plandelays(:,2)/C0*1e9;  % in ns
    if CORREL
        expdelays = TrigCor(ind,tag==1)';
    else
        expdelays = TrigTime(ind,tag==1)';  
    end
    expdelays = expdelays/FSAMPLING*1e9; 
    chi2ndfp(ind) = sum((expdelays-plandelays(:,2)).^2)./(errort.^2)/(multp(ind)-1);
    [LinParametersPlan Res] = FitSlope(plandelays(:,2),expdelays,1,0);
    slopep(ind) = LinParametersPlan(1);
    %
end;

if Type==0
    filename = [TEXT_PATH sprintf('R%d_sphrecons.txt', nrun)];
    disp(sprintf('Opening file %s...',filename)); 
    if fopen(filename)>0
        sphRes = load( filename );    
    else
        disp(sprintf('File %s does not exist.',filename));
        return
    end  
    disp 'Done.'
    
    for i=1:size(sphRes,1)
        
        ind = find(Struct.Coinc.IdCoinc==sphRes(i,1));
        
        flags(ind)=1;
        times(ind) = sphRes(i,2);
        mults(ind) = sphRes(i,3);
        x0(ind) = sphRes(i,4);
        y0(ind) = sphRes(i,5);
        z0(ind) = sphRes(i,6);
        if z0(ind)<REFALT
            z0(ind)=REFALT+(REFALT-z0(ind));
        end; 
        [rhos(ind) thetas(ind) phis(ind)]=Convert2Sph(x0(ind),y0(ind),z0(ind));
        %
        Xs=[x0(ind) y0(ind) z0(ind)];  % reconstructed source position
        tag = Tag(ind,:);
        indsci=find([Struct.Setup.Det.isScint]==1);
        tag(indsci)=0;  % put sci tag at 0 for radio reconstruction
        detId=DetId(ind,tag==1);
        detPos = [X(tag==1)' Y(tag==1)' Z(tag==1)'];
        %size(sqrt( sum( ( detPos - ones( mults(ind), 1 )*Xs ).^2, 2 ) ))
        sphdelays = [ detId' sqrt( sum( ( detPos - ones( mults(ind), 1 )*Xs ).^2, 2 ) )];  % in meters
        sphdelays(:,2) = sphdelays(:,2) - min(sphdelays(:,2));
        sphdelays(:,2) = sphdelays(:,2)/C0*1e9;  % in ns
        if CORREL
            expdelays = TrigCor(ind,tag==1)';
        else
            expdelays = TrigTime(ind,tag==1)';  
        end
        expdelays = expdelays/FSAMPLING*1e9; 
        chi2ndfs(ind) = sum((expdelays-sphdelays(:,2)).^2)./(errort.^2)/(multp(ind)-1);
        [LinParametersPlan Res] = FitSlope(sphdelays(:,2),expdelays,1,0);
        slopes(ind) = LinParametersPlan(1);
    end
end

%% Load dst
% dstname = sprintf('Fdst%d',nrun);
% dstname=[DST_PATH dstname];
% dst = load(dstname);
%Struct = dst.Struct;

%% Write to structure
disp('Now writting reconstruction results to DST...')
planRecons.Flag = flag;
planRecons.L = multp;
planRecons.T = timep;
planRecons.Theta = thetap;
planRecons.dTheta = dThetap;
planRecons.Phi = phip;
planRecons.dPhi = dPhip;
planRecons.Chi2 = chi2;
planRecons.Signif = signif;
planRecons.Chi2Delay = chi2ndfp;
planRecons.SlopeDelay = slopep;

if Type==0
    if ~isfield(Struct.Coinc,'DelayCorrRecons')
        Struct.Coinc.PlanRecons.Radio = planRecons;
    else
        Struct.Coinc.DelayCorrRecons.PlanRecons.Radio=planRecons;
    end;

    sphRecons.Flag = flag;
    sphRecons.L = mults;
    sphRecons.T = times;
    sphRecons.Rho = rhos;
    sphRecons.Theta = thetas;
    sphRecons.Phi = phis;
    sphRecons.X0 = x0;
    sphRecons.Y0 = y0;
    sphRecons.Z0 = z0;
    sphRecons.Chi2Delay = chi2ndfs;
    sphRecons.SlopeDelay = slopes;
    
    if ~isfield(Struct.Coinc,'DelayCorrRecons')
        Struct.Coinc.SphRecons = sphRecons;
    else
        Struct.Coinc.DelayCorrRecons.SphRecons=sphRecons;
    end;
    
elseif Type==1
    if ~isfield(Struct.Coinc,'DelayCorrRecons')
        Struct.Coinc.PlanRecons.Hybrid = planRecons;
    else
        Struct.Coinc.DelayCorrRecons.PlanRecons.Hybrid=planRecons;
    end;
end;

% %% Save dst
% filename = [DST_PATH sprintf(dst_filename,nrun,NbIterDst);
% save(filename,'Struct');
% disp 'Done.'
