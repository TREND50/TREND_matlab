function [delerror]=PlotDelays(nrun,cid,dst)
% Copied from soft/plotDelays.m
% Perform the delay plot for the reconstructed coincs
% 22/12/2010 OMH

% 
SharedGlobals;
dis = 1;

%% Load dst
if ~exist('dst')
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    dst = load(dstname);
end
Struct = dst.Struct;
Struct = Dist2Source(Struct);
CoincStruct = Struct.Coinc;
SetupStruct = Struct.Setup;
DetectorType = [dst.Struct.Setup.Det.isScint];
Detectors=[SetupStruct.Det.Name];
X = [SetupStruct.Det.X];
Y = [SetupStruct.Det.Y];
Z = [SetupStruct.Det.Z];
%

CoincId = CoincStruct.IdCoinc;
if ~exist('cid')   
    cid = CoincId;
end

L = CoincStruct.Mult;
Lant = CoincStruct.MultAnt;
%Lsci = CoincStruct.MultSci;
thetap = CoincStruct.PlanRecons.Radio.Theta;
phip = CoincStruct.PlanRecons.Radio.Phi;
dthetap = CoincStruct.PlanRecons.Radio.dTheta;
dphip = CoincStruct.PlanRecons.Radio.dPhi;
chi2p = CoincStruct.PlanRecons.Radio.Chi2;
signifp = CoincStruct.PlanRecons.Radio.Signif;

x = CoincStruct.SphRecons.X0;
y = CoincStruct.SphRecons.Y0;
z = CoincStruct.SphRecons.Z0;
rhos = CoincStruct.SphRecons.Rho;
thetas = CoincStruct.SphRecons.Theta;
phis = CoincStruct.SphRecons.Phi;
dist = Struct.Coinc.SphRecons.minDistSource;
%delerror = zeros(length(cid),Detectors);

%% Loop on coincs
for i=1:length(cid)  % Scan all coincs in the list
  id = cid(i); % Coinc number
  ind = find(CoincId==id);
  if isempty(ind)
      disp(sprintf('Coinc %d not found in list of reconstructed coincs.',id))
      continue
  end
  in = find(CoincStruct.Det.Tag(ind,:)==1);
  
  
  %% Experimental delays
  expdelays = CoincStruct.Det.TrigTime(ind,in)'/FSAMPLING*1e9;  % Measured delays (samples->ns)
  
  %% Compute expected delays
  %detId = CoincStruct.Det.Id(ind,in)';
  %[a ind_det]=intersect(Detectors,detId);
  ind_det = in
  detId = Detectors(in)'
  
  %if sum(DetectorType(ind_det))>0  % Skip scintillators
  %    continue
  %end
  xPos = X(ind_det)';
  yPos = Y(ind_det)';
  zPos = Z(ind_det)';
  detPos = [xPos yPos zPos];

  % Plane recons
  k= [sind(phip(ind))*sind(thetap(ind)) -cosd(phip(ind))*sind(thetap(ind)) -cosd(thetap(ind))]';  % Warning... With our X-Y convention, cartesian vector coordinates are different from usual 
  plandelays = [ detId detPos*k ];  % In meters
  plandelays(:,2) = plandelays(:,2) - min(plandelays(:,2)); % Use 1st triggering antenna as ref
  plandelays(:,2) = plandelays(:,2)/C0*1e9;  % in ns
  [a isAnt b] = intersect(plandelays(:,1),Detectors(DetectorType==0));
  [a isSci b] = intersect(plandelays(:,1),Detectors(DetectorType==1));
  
  % Spherical recons
  Xs=[x(ind) y(ind) z(ind)];  % reconstructed source position
  sphdelays = [ detId sqrt( sum( ( detPos - ones( L(ind), 1 )*Xs ).^2, 2 ) )];  % in meters
  sphdelays(:,2) = sphdelays(:,2) - min(sphdelays(:,2));
  sphdelays(:,2) = sphdelays(:,2)/C0*1e9;  % in ns
  
  plandelays_ant = plandelays(isAnt,:);
  plandelays_sci = plandelays(isSci,:);
  sphdelays_ant = sphdelays(isAnt,:);
  sphdelays_sci = sphdelays(isSci,:);
  expdelays_ant = expdelays(isAnt,:);
  expdelays_sci = expdelays(isSci,:);
  
  %delerror(ind,ind_det) = (sphdelays_ant(:,2)-expdelays_ant)/5;
  
  %% Plot
  errort = 1*ErrorTrig/FSAMPLING*1e9;   %ns
  [LinParametersPlan Res] = FitLin(plandelays_ant(:,2),expdelays_ant,[1,0],0);
  slope_plan = LinParametersPlan(1);
  khi2 = Res/errort.^2;
  khi2n_plan = khi2/Lant(ind);
  [LinParametersSph Res] = FitLin(sphdelays_ant(:,2),expdelays_ant,[1,0],0);
  slope_sph = LinParametersSph(1);
  khi2 = sum((sphdelays_ant(:,2)-expdelays_ant).^2/errort.^2);
  khi2n_sph = khi2/Lant(ind);
  
  % To be done: write Chi2 to struct.
  
  if dis
    disp(sprintf('Plane reconstruction: theta = %3.1f pm %3.1f deg, phi = %3.1f pm %3.1f deg, Chi2 = %3.1f, Signif = %3.1f',thetap(ind),dthetap(ind),phip(ind),dphip(ind),chi2p(ind),signifp(ind)));
    disp(sprintf('Spherical reconstruction: theta = %3.1f deg, phi = %3.1f deg, R = %5.1f m',thetas(ind),phis(ind),dist(ind)));
    disp(sprintf('Spherical reconstruction: x = %3.1f m, y = %3.1f m, z = %5.1f m',x(ind),y(ind),z(ind)));

    error = ones(Lant(ind),1)*errort;
    figure(2)
    title = sprintf('Delay plots Coinc %d Run %d',id,nrun);
    set(2,'Name',title,'NumberTitle','off')
    subplot(1,2,1)
    hold off
    errorbar(plandelays_ant(:,2), expdelays_ant, error ,'sk','MarkerFace','k' )
    hold on
    plot(plandelays_sci(:,2), expdelays_sci, 'sr','MarkerFace','r' )
    text(plandelays(:,2) + (max(plandelays(:,2))/15),expdelays,num2str(plandelays(:,1)),'FontSize',14);
    grid on
    %title(sprintf('Plane analysis : theta = %4.1f deg, phi = %4.1f deg',thetap(id),phip(id)));
    xlabel('Expected Delay [ns]',labelOpts{:} )
    ylabel('Measured Delay [ns]',labelOpts{:} )
    maxi = max(max(plandelays_ant(:,2)),max(expdelays_ant))*1.1-100;
    line([-100 maxi],[-100 maxi])
    text(-70,0.95*maxi,'Plane wave reconstruction','FontSize',12,'FontWeight', 'bold');
    text(-70,0.9*maxi,sprintf('theta=%3.1f deg, phi=%3.1f deg',thetap(ind),phip(ind)));
    text(-70,0.85*maxi,sprintf('Slope = %3.2f',slope_plan));
    text(-70,0.8*maxi,sprintf('chi^2/ndf = %3.2f',khi2n_plan));
    axis([-100 maxi -100 maxi]);
    %
    subplot(1,2,2);
    hold off
    errorbar(sphdelays_ant(:,2), expdelays_ant, error ,'sk','MarkerFace','k' )
    hold on
    plot(sphdelays_sci(:,2), expdelays_sci, 'sr','MarkerFace','r' )
    text(sphdelays(:,2) + (max(sphdelays(:,2))/10),expdelays,num2str(sphdelays(:,1)),'FontSize',14);
    grid on
    %title(sprintf('Spherical analysis : rho = %4.1f m, theta= %4.1f deg, phi = %4.1f deg',rhos(id),thetas(id),phis(id)));
    xlabel('Expected Delay [ns]',labelOpts{:} )
    ylabel('Measured Delay [ns]',labelOpts{:} )
    maxi = max(max(sphdelays(:,2)),max(expdelays))*1.1;
    line([-100 maxi],[-100 maxi])
    text(-70,0.95*maxi,'Spherical wave reconstruction','FontSize', 12,'FontWeight', 'bold');
    text(-70,0.9*maxi,sprintf('R=%3.1f m',dist(ind)));
    text(-70,0.85*maxi,sprintf('theta=%3.1f deg, phi=%3.1f deg',thetas(ind),phis(ind)));
    text(-70,0.8*maxi,sprintf('Slope = %3.2f',slope_sph));
    text(-70,0.75*maxi,sprintf('chi^2/ndf = %3.2f',khi2n_sph));  
    axis([-100 maxi -100 maxi]);
    
    %PlotCoinc(nrun,cid)
    if length(cid)>1
        pause
    end
  end
end


