function [] = generateDummyEvent(rhos,thetas,phis,id)
% Generate timetable for given source position (rhos,thetas,phis)
% to test reconstruction
% OMH 14/02/2014

SharedGlobals;
random = 1;
mult = 6;
%detsIn = [101:138 140];
AllDets = [101:110];
detsIn = zeros(1,mult);
for i = 1:mult
    ind = ceil(length(AllDets)*rand);
    detsIn(i) = AllDets(ind);
    AllDets(ind) = [];
end
detsIn =sort(detsIn);

%% Define source
[x0 y0 z0]=Convert2Cart(rhos,thetas,phis);
Xs=[x0 y0 z0];  % reconstructed source position

%% Load antenna positions
filename= 'coord_antennas_all.txt';
if fopen(filename)<0
    disp 'No antenna position file. Abort.'
    return
end
detPos = load(filename);
detId = detPos(:,1);
X = detPos(:,2);
Y = detPos(:,3);
Z = detPos(:,4);
[c indd] = intersect(detId,detsIn);

mult = length(indd);
X = X(indd);
Y = Y(indd);
Z = Z(indd);
detPos = [X Y Z];

%% Compute trigger times
delays = sqrt( sum( ( detPos - ones( mult, 1 )*Xs ).^2, 2 ) );  % in meters
delays = delays/C0*FSAMPLING;  % in sample units
err = ErrorTrig*randn(mult,1); %Error in samples
if random
    delays = delays + err;
end
delays = delays - min(delays);
[delays ind] = sort(delays);
detIdsorted = detId(ind);
txtTable = [zeros(mult,1) detIdsorted ones(mult,1) ones(mult,1) delays delays zeros(mult,1) zeros(mult,1) zeros(mult,1) zeros(mult,1) zeros(mult,1) zeros(mult,1)];

%% Write to file
filename = [TEXT_PATH sprintf('R%d_coinctable.txt',id)];
fid = fopen( filename, 'w' );        
for l = 1:size( txtTable, 1 )
    fprintf( fid, '%20f ', txtTable( l, 1 ) );   % Event time after correction (in sample counts)
    fprintf( fid, '%6d ',  txtTable( l, 2) );    % Antenna ID 
    fprintf( fid, '%6d ',  txtTable( l, 3) );    % Event number
    fprintf( fid, '%6d ',  txtTable( l, 4) );    % Coinc number
    fprintf( fid, '%8.2f ', txtTable( l, 5 ) );  % Delay with 1rst antenna (standard method)
    fprintf( fid, '%8.2f ', txtTable( l, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
    fprintf( fid, '%6.2f ', txtTable( l, 7 ) );  % Average correlation for this antenna (intercorrelation method)
    fprintf( fid, '%6d ', txtTable( l, 8 ) );  %  Saturation flag
    fprintf( fid, '%6.2f ', txtTable( l, 9 ) );  % Amplitude
    fprintf( fid, '%6.3f ', 0 ); % Gain OBSOLETE
    fprintf( fid, '%6.3f ', 0);  % Amplitude method 1 (bline)
    fprintf( fid, '%6.3f ', 0);  % Amplitude method 2 (PSD)
    fprintf( fid, '\n' );
end
fclose(fid);
disp(sprintf('Coinctable written to %s',filename))

%% Build antenna coordinates list 
filename = [TEXT_PATH sprintf('coord_antennas_%d.txt',id)];
fid = fopen(filename, 'w+' );
for kk = min(detId):max(detId)
    j = find( txtTable(:,2) == kk );
    if isempty( j )
        fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
        kk, 0, 0, 0 );
    else
        j = find(detId == kk );
        fprintf( fid, '%3d %6.1f %6.1f %6.1f\n', ...
        kk, X(j), Y( j ), Z( j ));
    end
end
fclose( fid );
disp(sprintf('Antennas coordinates written to %s',filename))
