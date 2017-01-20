function [] = loopAnaLog(periodId)
% Loop on ANaLog for all antennas and all runs in periodId
% Computes live time & T0 + T1 infos

SharedGlobals;
Antennas = [101:138 140 148:158];

if periodId>0
    switch (periodId)
      case 1
        nruns = 2538:2585;
      case 2
        nruns = 2685:2890;
      case 3
        nruns = 3000:3086;
      case 4
        nruns = 3157:3256;        
      case 5
        nruns = 3336:3371;   
      case 6
        nruns = 3562:3733;
        %nruns = 3577
      case 7
 %       nruns = 3835:3999;
        nruns = 3941:3999;
      case 8
        nruns = 4183:4389;
      case 9
        nruns = 4444:5014;
      case 10
%        nruns = 5070:5913;
        nruns = 5071:5913;
    end
else
    nruns =nrun;  % Provide run list by hand
end

%% Write to file
filename = sprintf('AnaLog_Period%d.txt',periodId);
nT0t = zeros(1,length(nruns));
nT1t = zeros(1,length(nruns));
livet = zeros(1,length(nruns));
for i = 1:length(nruns)
    nrun = nruns(i);
    livethisRun = zeros(1,length(Antennas));
    for j = 1:length(Antennas)
        ant = Antennas(j);
        [r f nT0 nT1 tlive] = AnaLog(nrun,ant,0);
        nT0t(i) = nT0t(i) + nT0;  % Sum up all T0s for this run
        nT1t(i) = nT1t(i) + nT1;  % Sum up all T1s for this run
        livethisRun(j) = tlive;
    end
    livet(i) = median(livethisRun);
    %
    fid = fopen( filename, 'a+' );        
    fprintf( fid, '%6d ', nrun);   % Event time after correction (in sample counts)
    fprintf( fid, '%15d ',  nT0t( i ) );    % Antenna ID 
    fprintf( fid, '%15d ',  nT1t( i ) );    % Event number
    fprintf( fid, '%13.2f ',  livet( i) );    % Coinc number
    fprintf( fid, '\n' );
    fclose(fid);
end
    
    
