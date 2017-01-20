function ok = checkTraj(theta,phi,xCore,yCore)
% checks if traj (theta,phi,xCore,yCore) crosses a mountain (ok = 0) or not (ok = 1)
% Developped in the framework of TREND50 EVA simus
% OMH 24/05/2013


SharedGlobals;
DISPLAY = 0;
ok = 1;

%% Compute traj in cart coord
rho = 1:10:20000;
st = sind(theta);
ct = cosd(theta);
sp = sind(phi);
cp = cosd(phi);

zCore = getElevation(xCore,yCore);
x = rho*st*sp+xCore;
y = rho*st*cp+yCore;
z = rho*ct+zCore;

if DISPLAY
    disp(sprintf('Shower [%d, %d] deg. Core position: (%d,%d) m',theta,phi,xCore,yCore))
    close all
    figure(1)
    subplot(3,1,1)
    plot(x,y)
    grid on
    xlabel('EW [m]', labelOpts{:})
    ylabel('SN [m]', labelOpts{:})
    subplot(3,1,2)
    plot(x,z)
    grid on
    xlabel('EW [m]', labelOpts{:})
    ylabel('Alt [m]', labelOpts{:})
    subplot(3,1,3)
    plot(y,z)
    grid on
    xlabel('SN [m]', labelOpts{:})
    ylabel('Alt [m]', labelOpts{:})
end

%% Get ground elev along traj
%disp 'Loading topology map...'
[xm,ym,zm] = getMap('map.bin',[min(x)-1000,max(x)+1000,min(y)-1000,max(y)+1000]);
%disp 'Done.'

if DISPLAY
  figure(2)
  contourf(xm,ym,zm,5)
  hold on
  xlabel('Easting [m]', labelOpts{:})
  ylabel('Northing [m]', labelOpts{:})
  colorbar
  plot(xCore,yCore','hm','MarkerFaceColor','m','MarkerSize',8)
  plot(x,y,'m','LineWidth',2)
end

zg = zeros(1,length(x));
for i=1:length(z)
    [a k] = min(abs(xm-x(i)));
    [a l] = min(abs(ym-y(i)));  
    zg(i) = zm(l,k);
end

below = find(z(rho>500)<zg(rho>500));
if length(below)>3  % Some points of the traj are below surface
    ok = 0;
end

if DISPLAY
    disp(sprintf('%d points below horizon.',length(below)))
    figure(3)
    plot(rho/1e3,zg,'k')
    hold on
    plot(rho/1e3,z,'r')
    title('Ground altitude along shower trajectory', labelOpts{:})
    xlabel('Distance to shower core [km]', labelOpts{:})
    ylabel('Altitude [m asl]', labelOpts{:})
    pause
end

