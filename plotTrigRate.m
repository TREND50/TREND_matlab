function plotTrigRate(runid,dst)
% Simply loads & plot TrigRate info from dsts
% OMH 15/07/2017

% disp 'in process...'
% return

SharedGlobals;

if ~exist('dst')
  dstname = [DST_PATH sprintf(dst_filename,runid,1)];
  fid2 = fopen(dstname);

    if fid2<0
        disp(sprintf('dst %s not found.',dstname))
    else
        fclose(fid2);
        disp(sprintf('Loading dst %s...',dstname))
        dst = load(dstname);
        disp 'Done.'
    end
end

%
t = dst.Struct.Setup.InfosRun.TrigTime; %[s]
if size(t,1) == 0
    disp 'No trigger rate measurment! Abort.'
    return
end
tstart = t(1);
tstop = t(end);
durs = (tstop-tstart);  % [s]
dur = durs/60; %[mn]
durh = dur/60; %[h]        
info = dst.Struct.Setup.InfosRun;
if ~isfield(info,'MetaTotalCoinc')
    disp 'No field MetaTotalCoinc' 
    nt1 = 0;
    ntf = 0;
else
    nt1 = sum(info.MetaTotalCoinc);
    ntf = sum(info.TotalCoincFiltered2);
end
%
coincrate = info.GlobalCoincRateRaw;
nt = size(coincrate,1);
if nt>10000
    disp(sprintf('Run %d: %3.1f h and %d points! Problem here. Abort.',runid,durh,nt))
    return
end
step = round(dur/nt); % [mn]
if step<1
    disp(sprintf('Error in computing trigger rate step: step = round(%3.1f/%d)=%d. Abort.',dur,nt,step))
    return
end
steph = step*60;
disp(sprintf('Run %d: %3.1f h, %d points (%3.2f / mn)',runid,durh,nt,step))
tmn = 0:step:dur;
figure(1)
subplot(1,1,1)
plot(tmn,coincrate(1:length(tmn)),'k','LineWidth',2)
xlim([0 max(tmn)])
xlabel('Time(mn)')
ylabel('T1 rate (Hz)')
grid on

detrate = info.DetCoincRateRaw;
det0rate = info.TrigRate;
size(det0rate)
nevts = [dst.Struct.Setup.Det.Evt];
dets = [dst.Struct.Setup.Det.Name];    
ndets = length(dets);

figure(2)
%subplot(1,2,2)
for i =1:4:ndets
    plot(tmn,det0rate(1:length(tmn),i)+i*20)
    xlim([0 max(tmn)])
    text(200,i*20+20,sprintf('A%d',dets(i)))
    hold on
    grid on
end
xlabel('Time(mn)')
ylabel('T1 rate (Hz)')
