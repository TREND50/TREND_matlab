function [Struct]=CppFileBuilder(Struct,opt)
% Select coincidences: good signals and good detector types.
% 06/12/10 OMH

if ~exist('opt')
    opt=0;
end

if opt==1  % Scintillators
    thresh = 3; % Recons for 3 scintillators at least
else  % Radio or Hybrid
    thresh = 4;  %Recons for 4 detectors at least
end

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
Evt=[RunSetup.Det.Evt];

CoincStruct = Struct.Coinc;
tag = CoincStruct.Det.Tag;
AmpMax = CoincStruct.Det.AmpMax;
TrigTime = CoincStruct.Det.TrigTime;
TrigCor = CoincStruct.Det.TrigCor;
CoefCor = CoincStruct.Det.CoefCor;

Gain = CoincStruct.Det.Gain;

ind_val = find(CoincStruct.IsValid==1);

%% Dump coordinates    
DumpCoordTxt(Struct);

%% Loop to selec coincs to be written to file
txtTable = [];
n = 0;
for j = 1:length(CoincStruct.IdCoinc)  % loop on all valid coincs
    %ind = ind_val(j);  % Index of the coincidence in the structure
    if opt ==2
        iDetsIn = find(tag(j,:)==1);  % All triggering dets with valid signals
    else
        iDetsIn =  find(tag(j,:)==1 & DetectorType==opt);  % Triggering antennas OR scints with valid signals
    end
    nDetsIn = size(iDetsIn,2);
    if nDetsIn>=thresh % radioRecons
        n = n+1;

        uTime = CoincStruct.Det.Time(j,iDetsIn)';
        CoincId = CoincStruct.IdCoinc(j)*ones(nDetsIn,1);
        AntId = CoincStruct.Det.Id(j,iDetsIn)';
        EvtId = CoincStruct.Det.Evt(j,iDetsIn)';
        trigTime = TrigTime(j,iDetsIn)';
        if CORREL
            correlTime = TrigCor(j,iDetsIn)';
            correlCoef = CoefCor(j,iDetsIn)';
        else
            correlTime = zeros(nDetsIn,1);
            correlCoef = zeros(nDetsIn,1);
        end
        flag = zeros(nDetsIn,1);
        Amp = AmpMax(j,iDetsIn)';
        gain = Gain(j,iDetsIn)';
                
        txtTable(end+1:end+nDetsIn,:) = [uTime AntId EvtId CoincId trigTime correlTime correlCoef flag Amp gain];
        
        decim = floor(length(ind_val)/10);
        if floor(j/decim)==j/decim
            disp(sprintf('Done at %2.0f percent',j/decim*10));
        end
    end
end
if size(txtTable,1)==0
    if opt==0
        disp(sprintf('No antenna coincs in R%d!',nrun));
    elseif opt==1
        disp(sprintf('No scintillator coincs in R%d!',nrun));
    elseif opt==2
        disp(sprintf('No hybrid coincs in R%d...',nrun));
    end
    return
end
txtTable = sortrows(txtTable,1);

%% Output to file
disp(sprintf('%d coincs to be written to R%d_coinctable.txt file.',n,nrun));
filename = [TEXT_PATH sprintf( 'R%d_coinctable.txt', nrun)];
fid = fopen( filename, 'w' );
for k = 1:size( txtTable, 1 )
  fprintf( fid, '%20f ', txtTable( k, 1 ) );   % Event time after correction (in sample counts)
  fprintf( fid, '%6d ',  txtTable( k, 2) );    % Antenna ID 
  fprintf( fid, '%6d ',  txtTable( k, 3) );    % Event number
  fprintf( fid, '%6d ',  txtTable( k, 4) );    % Coinc number
  fprintf( fid, '%8.2f ', txtTable( k, 5 ) );  % Delay with 1rst antenna (standard method)
  fprintf( fid, '%8.2f ', txtTable( k, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
  fprintf( fid, '%6.2f ', txtTable( k, 7 ) );  % Average correlation for this antenna (intercorrelation method)
  fprintf( fid, '%6d ', txtTable( k, 8 ) );  %  Noisy period flag
  fprintf( fid, '%6d ', txtTable( k, 9 ) );  % Amplitude
  fprintf( fid, '%6.3f ', txtTable( k, 10 ) );  % Gain
  fprintf( fid, '\n' );
end
fclose( fid );
clear txtTable
