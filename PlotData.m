function [  ] = PlotData( nrun,det,evt )
% Quickly plot data for detector det and run nrun
% OMH 15/02/2012

SharedGlobals;
DISPLAY = 1;

fd = OpenFileData( nrun, det);
if fd==-1
    disp(sprintf('Could not find data file for detector %d.',det))
    return
end

ft = OpenFileTime( nrun, det);
time=double(fread(ft,inf,'uint32'));

t = [0:ibuff-1]*5e-9*1e6;
maxi = zeros(1,length(evt));
sig = zeros(1,length(evt));
count_prev = 0;

%% Loop
for i = 1 :length(evt)
  fseek(fd,ibuff*(evt(i)-1),'bof');  % Set pointer at beginning of data
  DataEvt=double(fread(fd,ibuff,'uint8'));
  
  ind_time=(evt(i)-1)*ibufft;
  unixs = time(ind_time+1);
  buf = time(ind_time+2);
  subbuf = time(ind_time+3);
  index = time(ind_time+4);
  count =(buf-1)*2.^28 + (subbuf-1)*ibuff + index;
  disp(sprintf('Run % d - Detector %d - Evt %d: \ntime = %3.0f samples (%d - %d - %d)',nrun,det,evt(i),count,buf,subbuf,index))
  dt(i) = (count-count_prev)/FSAMPLING;
  disp(sprintf('Delay with previous event: %3.6f s',dt(i)))
  
  count_prev = count;
  m = mean(DataEvt);
  s = std(DataEvt);
  noise = find(abs(DataEvt-m)<s);
  maxi(i) = max(abs(DataEvt-m));
  stdnoise = std(DataEvt(noise));
  sig(i) = stdnoise;
  
  %%
  if DISPLAY 
      disp(sprintf('Detector %d - std dev noise = %3.3f ',det,stdnoise))
      figure(1)
      set(1,'NumberTitle','off','Name',sprintf('Detector %d - Event %d',det,evt(i)))
      subplot(2,1,1)
      plot(t,DataEvt,'k','Linewidth',1)
      grid on
      hold on
      plot(t,DataEvt,'k','Linewidth',1)
      plot(t(index),DataEvt(index),'ro','Linewidth',2)
      xlim([0 max(t)])
      xlabel('Time [mus]',labelOpts{:})
      ylabel('Amplitude [LSB]',labelOpts{:})
      hold off
      subplot(2,1,2)
      hist(DataEvt,100)
      pause
  end
end

figure(det)
subplot(2,1,1)
hist(maxi,100)
xlabel('Max amplitude above bline [LSB]', labelOpts{:})
subplot(2,1,2)
hist(dt,500)
xlabel('Delay with previous trigger [s]', labelOpts{:})


