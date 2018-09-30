function [] = loopScanCandidates(periodID)

[thpel phpel mult radel rampel run coinc val t] = ScanCandidateDST(periodID,0);
close all;
sel = find(mult>8 & (thpel>300 | thpel<60));
[run(sel)' coinc(sel)']

mult(sel)
for i = 2:length(run(sel))
    AnaCoinc(run(sel(i)),coinc(sel(i)))
    pause
    close all
end
return

pause
t = (t-t(1))/3600/24;
dt = diff(t);
dth = dt*24;
figure(12)
subplot(2,1,1)
plot(t,thpel,'+','MarkerSize',4)
grid on
xlabel('Time (day)')
ylabel('Zenith (deg)')
subplot(2,1,2)
plot(t,phpel,'+','MarkerSize',4)
grid on
xlabel('Time (day)')
ylabel('Azimuth (deg)')
% 
% figure(13) 
% subplot(2,1,1)
% hist(dt,100)
% datastats(dt)
% xlabel('\Deltat with previous candidate (days)')
% subplot(2,1,2)
% hist(dth(dth<1),100)
% datastats(dth(dth<1))
% xlabel('\Deltat with previous candidate (hours)')

pause
