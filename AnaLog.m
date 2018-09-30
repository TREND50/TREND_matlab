function [r f nT0 nT1 tlive] = AnaLog( nrun, antenna, dis )
% Read log file for antenna ant in nrun ant
% After code cleaning by Valentin
% Returns average trigger Rate over run & Fraction of events/spikes
% OMH 15/02/2012

if nrun<3561
    disp(sprintf('Error: log files format invalid before R3561... Launching AnaLogOldFormat(%d)',nrun))
    pause(3)
    AnaLogOldFormat(nrun,antenna);
    return
end

if ~exist('dis')
    dis = 1;
end
SharedGlobals;
j = 0;
r = 0;
f = 0;
nT0 = 0;
nT1 = 0;
tlive = 0;
for i = 1:length(antenna)
filename = sprintf( 'R%06d_A%04d_log.txt', nrun, antenna(i) ); 
filename = [LOG_PATH, filename];
if fopen(filename)<0
    disp(sprintf('Could not find file %s',filename))
    continue
    %return
end
logfile=load(filename);
format long
if size(logfile,1)==0
	disp 'Empty log file! Abort.'	
	continue
end

time_step = 2^28*5e-9; % time to write a buffer
tloop = logfile(:,1); % Time at begining of reading loop
tidle =  logfile(:,2); % Waiting time after reading 1/2 buffer 
tanat0 =  logfile(:,3); % Time to scan for spikes (T0) 
tanat1 =  logfile(:,4); % Time to scan for coincs (T1)
twrite =  logfile(:,5); % Time to write data to file
tprocess = tanat0 + tanat1 + twrite;
ttot = tidle + tprocess;
iloop = logfile(:,6);
irqs =  logfile(:,7);
irqe =  logfile(:,8);
nspikeb = logfile(:,9);
ntrigb = logfile(:,10);
lsb = logfile(:,11);
nspiket = cumsum(nspikeb);
ntrigt = cumsum(ntrigb);
if nspiket(end) == 0
    disp 'No spikes for this detector.'
    %return 
end
tmn = (tloop-tloop(1))/60;
spikerate = nspikeb/time_step;
trigrate = ntrigb/time_step;

check = irqe-irqs;
nlate = length(find(check>=1));
nverylate = length(find(check>2));
durlate = time_step*sum(nlate)/60; % minutes
    
%% Display
if dis
    figure(1)
    set(1,'NumberTitle','off','Name','Std dev evolution')
    plot(tmn,lsb,'k','Linewidth',2)
    xlabel('Time [mn]',labelOpts{:})
    ylabel('Std dev [LSB]',labelOpts{:})
    xlim([0 max(tmn)])
    grid on
    
    figure(2)
    set(2,'NumberTitle','off','Name','T0 rate')
%    subplot(2,1,1)
    plot(tmn,spikerate+j*200,'Linewidth',1)
    xlabel('Time (mn)',labelOpts{:})
    ylabel('T0 rate (Hz)',labelOpts{:})
    xlim([0 max(tmn)])
    text(100,j*200+100,sprintf('A%d',antenna(i)))
    grid on
    hold on
%     subplot(2,1,2)
    figure(3)
    set(3,'NumberTitle','off','Name','T1 rate')
%    subplot(2,1,1)
    plot(tmn,trigrate+j*200,'Linewidth',1)
    xlabel('Time (mn)',labelOpts{:})
    ylabel('T1 rate (Hz)',labelOpts{:})
    xlim([0 max(tmn)])
    text(100,j*200+100,sprintf('A%d',antenna(i)))
    grid on
    hold on

    %     plot(tmn,trigrate,'g','Linewidth',2)
%     xlabel('Time [mn]',labelOpts{:})
%     ylabel('T1 trig rate [Hz]',labelOpts{:})
%     grid on
%     xlim([0 max(tmn)])

    figure(30)
    set(30,'NumberTitle','off','Name','Event number')
    %subplot(2,1,1)
    plot(tmn,nspiket,'k','Linewidth',2)
    hold on
    plot(tmn,ntrigt,'g','Linewidth',2)
    xlabel('Time [mn]',labelOpts{:})
    ylabel('Number of spikes',labelOpts{:})
    xlim([0 max(tmn)])
    grid on
    hold off
    % subplot(2,1,2)
    % plot(tmn,ntrigt,'g','Linewidth',2)
    % xlabel('Time [mn]',labelOpts{:})
    % ylabel('Number of events',labelOpts{:})
    % grid on
    % xlim([0 max(tmn)])
    %return
    
    if 0
    figure(6)
    set(6,'Name','Wait time')
    subplot(2,1,1)
    plot(tmn,tidle,'+k')
    subplot(2,1,2)
    hist(tidle,100)
    figure(7)
    set(7,'Name','T0 scan time')
    subplot(2,1,1)
    plot(tmn,tanat0,'+k')
    subplot(2,1,2)
    hist(tanat0,100)    
    figure(8)
    set(8,'Name','T1 scan time')
    subplot(2,1,1)
    plot(tmn,tanat1,'+k')
    subplot(2,1,2)
    hist(tanat1,100)    
    figure(9)
    set(9,'Name','data write time')
    subplot(2,1,1)
    plot(tmn,twrite,'+k')
    subplot(2,1,2)
    hist(twrite,1000)    
    figure(10)
    set(10,'Name','Process time')
    subplot(2,1,1)
    plot(tmn,tprocess,'+k')
    subplot(2,1,2)
    hist(tprocess,100)    
    figure(11)
    set(11,'Name','Loop time')
    subplot(2,1,1)
    plot(tmn,ttot,'+k')
    line([0 max(tmn)],[time_step time_step],'Color','r')
    subplot(2,1,2)
    hist(ttot,100)  
    figure(12)
    set(12,'Name','Check')
    subplot(2,1,1)
    plot(tmn,check,'+k')
    subplot(2,1,2)
    hist(check,100)  
    end
    j = j+1
end

nT0 = nspiket(end);
nT1 = ntrigt(end);
r = nT1/tmn(end)/60;
f = nT1/nT0*100;
tlive = tmn(end)-durlate;
disp(sprintf('Run %d Antenna %d - %3.1f mins', nrun, antenna(i), tmn(end)))
disp(sprintf('%d spikes (av rate = %3.1f Hz) & %d events recorded (av rate = %3.1f Hz). Ratio = %3.2f pc.',nT0,mean(spikerate),nT1,mean(trigrate),f))
disp(sprintf('Out of sync for %d buf(s) ( = %3.2f mn, %3.2f pc of acq time)',nlate,  durlate, durlate/tmn(end)*100))
fclose all;

end

