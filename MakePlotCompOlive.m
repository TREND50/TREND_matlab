function [] = MakePlotComp()

SharedGlobals;

TXTPATH = '../data/candidates/candidates_dc';

E={'8e16','1e17','2e17','3e17','5e17','7e17','1e18','3e18'}
%wini = [ 28.2476   46.4490   39.2201   38.2236   31.5862   17.4193   13.0303    1.2840] % erf fit
%wini = [24.9828   39.3586   31.1448   29.4424   25.9159   16.2594   12.9954    1.2840];  %erf fit 10% E bias
wini = [23.4766   36.1732   27.6280   25.5481   23.0133   15.4282   12.9495    1.2840]; %erf fit 15% error
%wini = [22.0498   33.2111   24.4344   21.9977   20.1427   14.4137   12.8653  1.2840];  %erf fit 20% E bias
%wini = [13.0991   63.3256   71.0977   54.2467   30.8900   17.7835   13.5920  1.4494];  %lin fit
%wini = [4.0250   44.3261   56.9189   45.6331   26.8680   16.1121   13.4927  1.4275];  % lin fit 10% Ebias
%wini = [0.0629   27.3740   43.9641   36.6492   23.0950   13.7707   13.3786    1.4056];  % lin fit 20% Ebias
wini = [28.2476   46.4490   39.2201   38.2236   31.5862   17.4193   13.0303    1.2840];  % erf fit 
wini = [23.4766   36.1732   27.6280   25.5481   23.0133   15.4282   12.9495    1.2840];  % erf fit 15% Ebias


%% Get DataChallenge results
for i = 1:length(E)
  selfile=load(sprintf('%s/%s/log_candanalysisend.txt',TXTPATH,E{i}), 'r' );
  out = find(selfile(:,2)>75 & selfile(:,3)<=350 & selfile(:,3)>=250);  % reject trajectories below hiron (theta >75° at SOuth)
  %out = [];
  disp(sprintf('E=%seV: killing %d/%d events below horizon',E{i},length(out),size(selfile,1)))
  in = setxor([1:size(selfile,1)],out);
  w(i) = wini(i)*length(in)/size(selfile,1);
  % Load valid trajectories
  theta=selfile(in,2);
  phi=selfile(in,3);
  tag=selfile(in,33);
  % Compute zen & az hists for each energy slice
  [count_theta(i,:),center_theta(i,:)]=hist(theta(tag==1),[5:10:85]);
  [count_phi(i,:),center_phi(i,:)]=hist(phi(tag==1),[20:40:340]);
  % Apply weigth fpr each E
  cth(i,:) = count_theta(i,:)./sum(count_theta(i,:))*w(i);
  cph(i,:) = count_phi(i,:)./sum(count_phi(i,:))*w(i);
end
totevt = sum(w);
disp(sprintf('There are %3.1f events expected from data challenge (%3.1f before horizon cut)',totevt,sum(wini)))

count_theta = sum(cth,1); % Sum weigthed hists
count_theta = count_theta*totevt/sum(count_theta); % Scale to DataChallenge nb of events 
count_theta_error = sqrt(count_theta);
count_theta_error = count_theta_error*totevt/sum(count_theta);
%
count_phi = sum(cph,1);
count_phi = count_phi*totevt/sum(count_phi);
count_phi_error = sqrt(count_phi);
count_phi_error = count_phi_error*totevt/sum(count_phi);

%% Now load experimental data
selfile=load(sprintf('Candidates_all.txt'), 'r' );  % All periods
theta=selfile(:,1);
phi=mod(selfile(:,2)+90+phigeom,360);   % North is Magnetic North in Zhaires

%% Theta plot
figure(1)
[countdata,center]=hist(theta,[5:10:85]);
countdata = countdata/sum(countdata)*204;  % Rescale to the 204 events expected in period 6
countdata_error = sqrt(countdata)/sum(countdata)*204;
fig=errorbar(center,countdata,countdata_error,'k-*','linewidth',2);
hold on
grid on
errorbar(center,count_theta,count_theta_error,'r-*','linewidth',2)
xlim([0 90])
xlabel('Theta [deg]','fontsize',14)
ylabel('Count [10 deg^{-1}]','fontsize',14)
set(gca, 'FontSize', 15)
saveas(fig,sprintf('thetacomp.png'))

figure(2)
[countdata,center]=hist(phi,[20:40:340]);
countdata = countdata/sum(countdata)*204;
countdata_error = sqrt(countdata)/sum(countdata)*204;
fig2 = errorbar(center,countdata,countdata_error,'k-*','linewidth',2);
hold on
grid on
errorbar(center,count_phi,count_phi_error,'r-*','linewidth',2)
xlim([0 360])
xlabel('Phi [deg]','fontsize',14)
ylabel('Count [40 deg^{-1}]','fontsize',14)
set(gca, 'FontSize', 15)
saveas(fig2,sprintf('phicomp.png'))

