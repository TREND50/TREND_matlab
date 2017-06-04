function [r f] = AnaLogOldFormat( nrun, antenna, dis )
% Read log file for antenna ant in nrun ant
% Before code cleaning by Valentin (Feb 2012)
% Returns average trigger rate over run 
% OMH 27/02/2012

r = 0;
f = 0;
% if nrun>=3561
%     disp(sprintf('Error: AnaLogOldFormat valid priori to R3561 only... Launching AnaLog(%d)',nrun))
%     pause(3)
%     AnaLog(nrun,antenna);
%     return
% end

if ~exist('dis')
    dis = 1;
end
SharedGlobals;
r = 0;
f = 0;
filename = sprintf( 'R%06d_A%04d_log.txt', nrun, antenna ); 
filename = [LOG_PATH, filename];
if fopen(filename)<0
    disp(sprintf('Could not find file %s',filename))
    return
end
logfile=load(filename);
format long
time_step = 2^28*5e-9/2; % Time to read a HALF buffer;
ts = logfile(:,1); % Time at begining of reading loop
strig = find(ts>0);
if size(strig,1)==0
    disp 'No data recorded on this detector'.
    return;
end

t_firsttrig = strig(1)*time_step*2;
irq =  logfile(:,2); % Buffer index
parity = [diff(irq)' 0]';
t = t_firsttrig + irq.*time_step*2+parity.*time_step;
tmn = (t-t(1))/60;  % Time in minutes since run start
nspike =  logfile(:,3); % number of spikes in this half buffer
ntrig =  logfile(:,4); % number of events recorded on disc in this half buffer
%rate = logfile(:,5); % Spike rate [Hz]
srate = nspike/time_step;
trate = ntrig/time_step;
snspike = cumsum(nspike);
sntrig = cumsum(ntrig);

ns = snspike(end);
nt = sntrig(end);
r = nt/tmn(end)/60;
f = nt/ns*100;
disp(sprintf('Run %d Antenna %d - %3.1f mins', nrun, antenna, tmn(end)))
disp(sprintf('%d spikes (av rate = %3.1f Hz) & %d events recorded (av rate = %3.1f Hz). Ratio = %3.2f pc.',ns,mean(srate),nt,mean(trate),f))

%% Plots
if dis
    figure(1)
    set(1,'Name', 'Event count','NumberTitle','off')
    plot(tmn,snspike,'k','LineWidth',2)
    hold on
    plot(tmn,sntrig,'g','LineWidth',2)
    grid on
    xlabel('Time [mn]', labelOpts{:})
    ylabel('Event nb', labelOpts{:})
    xlim([0 tmn(end)])
    %
    figure(2)
    set(2,'Name', 'Trig rate','NumberTitle','off')
    subplot(2,1,1)
    plot(tmn,srate,'k','LineWidth',2)
    xlim([0 tmn(end)])
    grid on
    xlabel('Time [mn]', labelOpts{:})
    ylabel('Spike rate [Hz]', labelOpts{:})
    subplot(2,1,2)
    plot(tmn,trate,'g','LineWidth',2)
    xlim([0 tmn(end)])
    grid on
    xlabel('Time [mn]', labelOpts{:})
    ylabel('Trig rate [Hz]', labelOpts{:})
end