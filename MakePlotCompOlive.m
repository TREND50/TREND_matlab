function [] = MakePlotComp()

SharedGlobals;

TXTPATH = '../data/candidates/candidates_dc';

E={'7e16','8e16','1e17','2e17','3e17','5e17','7e17','1e18','2e18','3e18'}
%wini = [4.3805   11.6370   35.0830   44.7381   57.8327   43.9752   19.0686   11.3971    2.9258    0.4473] % erf fit
%wini = [3.1023    7.9324   21.3440   23.2440   30.2880   30.7937   18.0414 11.3921    2.9258    0.4473];  %erf fit (f=0.8)
wini = [0    1.4767   16.1715   31.1369   33.2288   23.2604   13.3545   10.8791    3.3457    0.5215]; %multivariate fit f=0.75
%wini = [0    3.0378   20.9780   37.7414   38.0718   25.5232   14.1998 11.2469   3.3716    0.5197] %multivariate fit f=0.8
%wini = [1.6181    6.9699   40.1474   60.6394   52.9439   31.7979   16.3448   12.0796    3.4063    0.5122] %multivariate fit f=0.95
%wini = [2.3371   11.9533   63.3876   70.6742   54.0201   30.6391   18.5805   11.7129    3.3851    0.5076];  %lin fit
%wini = [0       0   28.0629   43.7333   36.5441   22.8290   13.9559   11.6511    3.2201    0.5312];  % lin fit 20% Ebias


%% Get DataChallenge results
for i = 1:length(E)
  selfile=load(sprintf('%s/%s/log_candanalysisend.txt',TXTPATH,E{i}), 'r' );
  tag = find(selfile(:,end)==1);  % Valid showers
  selfile = selfile(tag,:);
  out = find(selfile(:,2)>70 & selfile(:,3)<=350 & selfile(:,3)>=250);  % reject trajectories below hiron (theta >75° at SOuth)
  %out = [];
  disp(sprintf('E=%seV: killing %d/%d events below horizon',E{i},length(out),size(tag,1)))
  in = setxor([1:length(tag)],out);
  w(i) = wini(i)*length(in)/length(tag);
  % Load valid trajectories
  theta=selfile(in,2);
  phi=mod(selfile(in,3)+phigeom-90,360);
  %tag=selfile(in,33);
  % Compute zen & az hists for each energy slice
  [count_theta(i,:),center_theta(i,:)]=hist(theta,[5:10:85]);
  [count_phi(i,:),center_phi(i,:)]=hist(phi,[20:40:340]);
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
%selfile=load(sprintf('ab.txt'), 'r' );  % All periods   Mult 5 & soft cuts
%selfile=load(sprintf('ab2.txt'), 'r' );  % All periods  Mult 6 & hard cuts
%selfile=load(sprintf('ab3.txt'), 'r' );  % All periods  Mult 5 & hard cuts
%selfile=load(sprintf('Candidates_Period6.txt'), 'r' );  % Period 6
theta=selfile(:,1);
%phi=mod(selfile(:,2)+90-phigeom,360);   % North is Magnetic North in Zhaires
phi=mod(selfile(:,2)-phigeom,360);   % North is Magnetic North in Zhaires
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

%% Phi plot
figure(2)
[countdata,center]=hist(phi,[20:40:340]);
countdata = countdata/sum(countdata)*204;
countdata_error = sqrt(countdata)/sum(countdata)*204;
fig2 = errorbar(center,countdata,countdata_error,'k-*','linewidth',2);
hold on
grid on
errorbar(center,count_phi,count_phi_error,'r-*','linewidth',2)
xlim([0 360])
set(gca,'xtick',[0:40:360])
xlabel('Phi [deg]','fontsize',14)
ylabel('Count [40 deg^{-1}]','fontsize',14)
set(gca, 'FontSize', 15)
saveas(fig2,sprintf('phicomp.png'))

