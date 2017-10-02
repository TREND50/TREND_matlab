function [] = DisplayTrigRate(runbeg,runend)
% Displays results from trigrate.mat
% Scipt used for presentation at 2015 GRAND workshop

SharedGlobals;
DISPLAY = 1;
ants = [101:138 140 148:158];
if ~exist('runend')
     runbeg = 2500;
     runend = 5070;   
end

%% Load structure
filename = ([MONITOR_PATH 'trigrate.mat']);
fid = fopen(filename);
if fid<0
    disp(sprintf('Could not fond %s. Abort',filename))
    return
else
    disp(sprintf('Loading %s',filename))
    fclose(fid);
    t = load(filename);
end

allruns = t.allruns;
nruns = size(allruns,2);
runs = [];
durh = [];
nt1 = [];
ntf = [];
detnball = zeros(1,length(ants));
detnbt1all = zeros(1,length(ants));
isin = zeros(1,length(ants));
for i =1:nruns
  if allruns{i}.id<runbeg
      continue
  elseif allruns{i}.id>runend
      break
  end
  runs(end+1) = allruns{i}.id; 
  durh(end+1) = allruns{i}.dur;
  nt1(end+1) = allruns{i}.nt1;
  ntf(end+1) = allruns{i}.ntf;
  %
  thisDets = allruns{i}.dets;
  det_coincrate = allruns{i}.detcrate;
  det_nbt1 = allruns{i}.detnbt1;
  [c ia ib] = intersect(ants,thisDets);
  thisDetRate = det_coincrate(ib);
  thisDetRate(~isfinite(thisDetRate)) = 0;
  %
  detnball(end+1,ia) = round(thisDetRate*durh(end)*3600); %Nb events on each antenna for this run (computed from T1 rate integration)
  detnbt1all(end+1,ia) = det_nbt1(ib); %Nb events on each antenna for this run (true value)
  isin(end+1,ia) = thisDetRate>0;
end    
detnball(1,:) = [];
isin(1,:) = [];
detnbt1all(1,:) = [];
eastclose = find(ants<=108 | ants==110 | ants==111);
ratet1_eastclose = sum(detnbt1all(:,eastclose).*isin(:,eastclose),2)./sum((isin(:,eastclose)>0),2)./durh'/3600;
rate_eastclose = sum(detnball(:,eastclose).*isin(:,eastclose),2)./sum((isin(:,eastclose)>0),2)./durh'/3600;
eastfar = find(ants>=112 | ants==109);
ratet1_eastfar = sum(detnbt1all(:,eastfar).*isin(:,eastfar),2)./sum((isin(:,eastfar)>0),2)./durh'/3600;
rate_eastfar = sum(detnball(:,eastfar).*isin(:,eastfar),2)./sum((isin(:,eastfar)>0),2)./durh'/3600;
west = find(ants>120 & ants<=132);
ratet1_west = sum(detnbt1all(:,west).*isin(:,west),2)./sum((isin(:,west)>0),2)./durh'/3600;
rate_west = sum(detnball(:,west).*isin(:,west),2)./sum((isin(:,west)>0),2)./durh'/3600;
cross = find(ants>=133);
ratet1_cross = sum(detnbt1all(:,cross).*isin(:,cross),2)./sum((isin(:,cross)>0),2)./durh'/3600;
rate_cross = sum(detnball(:,cross).*isin(:,cross),2)./sum((isin(:,cross)>0),2)./durh'/3600;
detnbt1sum = sum(detnbt1all(runs>=runbeg & runs<=runend,:),1);

%% Plots
[unixs iasc] = sort(t.unixs);
[yi, mi, di, hi, mni, si] = UnixSecs2Date(unixs);
d = datenum(yi,mi,di,hi,mni,si);
t0rate = t.t0rate(:,iasc);

ind = 0;
for i=1:length(ALLDETS)
    poss = floor(i/20.1);
    if poss~=floor((i-1)/20.1)
        ind = 0;
    end
    ind=ind+1;
    figure(poss+1)
    subplot(4,5,ind)
    plot(d,t0rate(i,:),'k-','MarkerSize',2)
    grid on
    datetick('x','mm/yyyy')
    xlabel(sprintf('Antenna %d',ALLDETS(i)))
    %xlabel('Date', labelOpts{:})
    %ylabel('Raw coinc rate [Hz]', labelOpts{:})
    sel = find(t0rate(i,:)>0);
    trun = t0rate(i,sel);
    nh = length(sel);
    disp(sprintf('Antenna %d: Fraction of time <10Hz: %3.2f pc, above 100Hz: %3.2f pc, above 500Hz: %3.2f pc', ALLDETS(i),length(find(trun<10))/nh*100,length(find(trun>100))/nh*100,length(find(trun>500))/nh*100))
    
end
pause
coincrate = t.coincrate(iasc);
smoothrate = smooth(coincrate,24*7);
%     unixs = t.unixs;
%     coincrate = t.coincrate;
figure(27)
in = find(detnbt1sum>0);
plot(ants(in),detnbt1sum(in),'sk','MarkerFaceColor','k')
grid on
xlabel('Antenna ID')
ylabel('Total T1s')
pause

figure(10)
plot(d,coincrate,'k')
hold on
plot(d,smoothrate,'g','LineWidth',2)
grid on
datetick('x','mm/yyyy')
xlabel('Date', labelOpts{:})
ylabel('Raw coinc rate [Hz]', labelOpts{:})
%
%Antenna group trig rate
ymax = max([max(ratet1_eastclose),max(ratet1_eastfar),max(ratet1_west),max(ratet1_cross)])*1.1;
figure(20)
subplot(2,2,1)
plot(runs,ratet1_eastclose,'k')
hold on
grid on
plot(runs,rate_eastclose,'k--')
ylim([0 ymax]);
xlabel('Run ID', labelOpts{:})
ylabel('CloseEast mean T1 rate [Hz]', labelOpts{:})
subplot(2,2,2)
plot(runs,ratet1_eastfar,'g')
hold on
grid on
plot(runs,rate_eastfar,'g--')
ylim([0 ymax]);
xlabel('Run ID', labelOpts{:})
ylabel('FarEast mean T1 rate [Hz]', labelOpts{:})
subplot(2,2,3)
plot(runs,ratet1_west,'b')
hold on
grid on
plot(runs,rate_west,'b--')
xlabel('Run ID', labelOpts{:})
ylabel('West mean T1 rate [Hz]', labelOpts{:})
ylim([0 ymax]);
subplot(2,2,4)
plot(runs,ratet1_cross,'r')
hold on
grid on
plot(runs,rate_cross,'r--')
ylim([0 ymax]);
xlabel('Run ID', labelOpts{:})
ylabel('Cross mean T1 rate [Hz]', labelOpts{:})
%
rfilt = ntf./nt1;
rf = ntf./durh/3600;
rt1 = nt1./durh/3600;
figure(23)
plot(runs,rt1,'k+','MarkerSize',2)
avrew = mean(rt1(runs<3920))
avrns = mean(rt1(runs>=3920))
sum(durh(runs<4389))/24;
sum(durh(runs>=4444))/24;
h = line([min(runs) 3920],[avrew avrew]);
set(h,'Color','r','LineWidth',4)
h = line([3920 max(runs)],[avrns avrns]);
set(h,'Color','r','LineWidth',4)

smoothrate = smooth(rt1,30);
hold on
plot(runs,smoothrate,'g','LineWidth',2)
ylim([0 max(rt1*1.1)])
grid on
xlabel('Run ID', labelOpts{:})
ylabel('Raw coinc rate [Hz]', labelOpts{:})
pause

figure(30)
subplot(3,1,1)
semilogy(runs,rfilt,'k')
grid on
xlabel('Run ID', labelOpts{:})
ylabel('Filt/Raw', labelOpts{:})
subplot(3,1,2)
plot(runs,rt1,'k')
ylim([0 max(rt1*1.1)])
grid on
xlabel('Run ID', labelOpts{:})
ylabel('Raw event rate [Hz]', labelOpts{:})
subplot(3,1,3)
plot(runs,rf,'k')
ylim([0 max(rf*1.1)])
grid on
xlabel('Run ID', labelOpts{:})
ylabel('Filt event rate [Hz]', labelOpts{:})

%% Displays
disp(sprintf('Slice %d-%d',min(runs),max(runs)))
disp(sprintf('Nb of runs = %d',length(runs)))
disp(sprintf('Total duration = %3.1f days',sum(durh)/24))
disp(sprintf('Mean coinc rate (before filtering) = %3.1f Hz.',mean(coincrate)))
disp(sprintf('Array trig rate (before filtering) = %3.1f Hz.',mean(rt1)));
disp(sprintf('Fraction of events passing filter = %d/%d = %3.3f.',sum(ntf),sum(nt1),sum(ntf)./sum(nt1)))
disp(sprintf('Array trig rate (after filtering) = %3.1f Hz.',mean(rf)));
disp(sprintf('Mean T1 trig rate on Close East group = %3.2f Hz.',mean(ratet1_eastclose(isfinite(ratet1_eastclose)))));
disp(sprintf('Mean T1 trig rate on Far East group = %3.2f Hz.',mean(ratet1_eastfar(isfinite(ratet1_eastfar)))));
disp(sprintf('Mean T1 trig rate on West group = %3.2f Hz.',mean(ratet1_west(isfinite(ratet1_west)))));
disp(sprintf('Mean T1 trig rate on Cross group = %3.2f Hz.',mean(ratet1_cross(isfinite(ratet1_cross)))));
disp(sprintf('Mean trig rate (after offline coinc search) on Close East group = %3.2f Hz.',mean(rate_eastclose(isfinite(rate_eastclose)))));
disp(sprintf('Mean trig rate (after offline coinc search) on Far East group = %3.2f Hz.',mean(rate_eastfar(isfinite(rate_eastfar)))));
disp(sprintf('Mean trig rate (after offline coinc search) on West group = %3.2f Hz.',mean(rate_west(isfinite(rate_west)))));
disp(sprintf('Mean trig rate (after offline coinc search) on Cross group = %3.2f Hz.',mean(rate_cross(isfinite(rate_cross)))));

