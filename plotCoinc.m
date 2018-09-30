function []=plotCoinc(nrun,IdCoinc)
% Adapted from PlotEvent
% OMH 30/09/2018

SharedGlobals;


%% Load dst
dstname = [DST_PATH sprintf(dst_filename,nrun,1)];   % Standard dst
dst = load(dstname);
Struct = dst.Struct;
Ncoinc=find(Struct.Coinc.IdCoinc==IdCoinc);
Tag=Struct.Coinc.Det.Tag(Ncoinc,:);
Evt=Struct.Coinc.Det.Evt(Ncoinc,:);
Detectors=[Struct.Setup.Det.Name];

in = find(Evt>0);
dets = Detectors(in);
Evt = Evt(in);
mult = length(dets);
DataEvt = {};


%% Get calib info
cal = load("calibCandidates.txt");
runId = cal(:,1);
sel = find(runId==nrun);
if length(sel)>0
  antCalib = cal(sel,2);
  gain = cal(sel,3);  % To go back to V/V gain 
  display 'Gain PSD:'
  [antCalib gain]    
else
  display('No gain available for this coinc')
  gain = ones(length(in));
  antCalib = dets;
end

%% Now get data
display(sprintf('R%dC%d:',nrun,IdCoinc))
for i=1:mult
  display(sprintf('Antenna %d Evt %d',dets(i),Evt(i)))
  fd = OpenFileData( nrun, dets(i));
  if fd==-1
    disp(sprintf('Could not find data file for detector %d.',dets(i)))
    return
  end
  fseek(fd,ibuff*(Evt(i)-1),'bof');  % Set pointer at beginning of data
  d=double(fread(fd,ibuff,'uint8'));
  inc = find(antCalib==dets(i));
  DataEvt{i} = (d-mean(d))/gain(inc)*1e6;  % mV
end

for i = 1:mult
  amax(i) = max(abs(DataEvt{i}));
end
amax = max(amax);

%% Now plot
t = 1:length(DataEvt{1});
t = t/FSAMPLING*1e6;
figure(1)
for i=1:mult
    subplot(3,3,i)
    plot(t,DataEvt{i})
    xlim([0 max(t)])
    ylim([-amax amax])
    xlabel('time (µs)')
    ylabel('Voltage (µV)')
    title(sprintf('Antenna %d',dets(i)))
end
