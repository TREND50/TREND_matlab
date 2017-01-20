function [antid dist ] = computeShowerDistance(theta,phi,posCore,posAnt)
% Compute distance to shower for all antennas in layout
% Shower caras (theta, phi, xcore [x, y, z]) following TREND conventions
% posAnt = [id x y z]
% OMH 19/11/2015

SharedGlobals;
DISPLAY = 0;

%% Load antenna position
if ~exist('posAnt')
    posAnt = load('coord_antennas_TREND50.txt');
    disp 'Loading antenna positions...'
    %posAnt = load('coord_antennas_GRANDproto.txt');
end
antid = posAnt(:,1);
posDet = posAnt(:,2:4);

%% Compute COre position
if length(posCore) == 2
    posCore(3) = getElevation(posCore(1),posCore(2));
end
posCore = (ones(length(antid),1)*posCore);

%% Compute distance to Core
dX =  posDet - posCore;
cp = cosd(phi); 
sp = sind(phi);
ct = cosd(theta); 
st = sind(theta);
u  = [ -sp*st, cp*st, ct ];  %% Warning!!!! This is different from usual since here x=WE & y=SN
U  = ones( length(antid), 1 )*u;
dist = sqrt( sum( dX.^2, 2 ) - sum( dX.*U, 2 ).^2 ); %Distance from antennas to shower axis   

if DISPLAY == 1
    in = find(posDet(:,3)>0);
    plot(antid(in),dist(in),'+k','MarkerSize',6)
    xlabel('Antenna ID')
    ylabel('Distance to shower [m]')
    grid on
end