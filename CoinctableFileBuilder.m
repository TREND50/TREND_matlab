function [Struct]=CoinctableFileBuilder(Struct,Type)
% Build table for C++ analysis
% 06/12/10 OMH

% AnalysisType=0 : radio only (MultAnt>=4, potential sci tag removed)
% AnalysisType=1 : hybrid
% AnalysisType=2 : scintillators

thresh=4; % in all cases

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
Evt=[RunSetup.Det.Evt];

CoincStruct = Struct.Coinc;
tag = CoincStruct.Det.Tag;
multant = CoincStruct.MultAnt;
%multsci = CoincStruct.MultSci;
indsci=find([RunSetup.Det.isScint]==1);
indant=find([RunSetup.Det.isScint]==0);
MinRaw = CoincStruct.Det.MinRaw;
MaxRaw = CoincStruct.Det.MaxRaw;
AmpMax = CoincStruct.Det.AmpMax;
TrigTime = CoincStruct.Det.TrigTime;
TrigCor = CoincStruct.Det.TrigCor;
DelayCor=zeros(1,length(Detectors));


% DelayCor(Detectors==101) = +5;
% DelayCor(Detectors==105) = +5;
% DelayCor(Detectors==106) = +7;
% DelayCor(Detectors==108) = +12;
% DelayCor(Detectors==113) = -5;
% DelayCor(Detectors==128) = -11;
% DelayCor(Detectors==133) = -5;
% DelayCor(Detectors==137) = +11;
% DelayCor(Detectors==138) = +11;
% DelayCor(Detectors==148) = +3;
% DelayCor(Detectors==149) = -3;
% DelayCor(Detectors==158) = -3;

if isfield(Struct.Coinc,'DelayCorrRecons')
    for i=1:length(CoincStruct.IdCoinc)
        indtag=indant(tag(i,indant)==1);
        %RunSetup.Det.DelayCorr(indtag)
        %TrigTime(i,indtag)
        for j=1:length(indtag)
%          	TrigTime(i,indtag(j))=TrigTime(i,indtag(j))-DelayCor(indtag(j)); % Units = samples
%             TrigCor(i,indtag(j))=TrigCor(i,indtag(j))-DelayCor(indtag(j));
        	TrigTime(i,indtag(j))=TrigTime(i,indtag(j))-RunSetup.Det(indtag(j)).DelayCorr;
            TrigCor(i,indtag(j))=TrigCor(i,indtag(j))-RunSetup.Det(indtag(j)).DelayCorr;
        end;         
       TrigTime(i,indtag) = TrigTime(i,indtag)-min(TrigTime(i,indtag));
       TrigCor(i,indtag) = TrigCor(i,indtag)-min(TrigCor(i,indtag));
        clear indtag
    end;
end;

CoefCor = CoincStruct.Det.CoefCor;

FlagTable=0;

ind_val = find(tag==1);

%% Create text file 
filename = [TEXT_PATH sprintf( 'R%d_coinctable.txt', nrun)];
fid = fopen( filename, 'w' );

%% Dump coordinates    
DumpCoordTxt(Struct);

%% Loop to select coincs to be written to file
n = 0;
for j = 1:length(CoincStruct.IdCoinc)  % loop on all valid coincs
    txtTable = [];
    
    if Type==0
        iDetsIn = indant(tag(j,indant)==1);
    else
        if multant(j)>0 & multsci(j)>0
            iDetsIn = find(tag(j,:)==1);
        else
            iDetsIn=0;
        end;
    end;
    nDetsIn = size(iDetsIn,2);

    if nDetsIn>=thresh % radioRecons
        
        FlagTable=1;
        
        n = n+1;

        %uTime = CoincStruct.Det.Time(j,iDetsIn)';
        uTime = CoincStruct.Det.UnixTime(j,iDetsIn)';
        CoincId = CoincStruct.IdCoinc(j)*ones(nDetsIn,1);
        %AntId = CoincStruct.Det.Id(j,iDetsIn)';
        AntId = Detectors(iDetsIn)';
        EvtId = CoincStruct.Det.Evt(j,iDetsIn)';
        trigTime = TrigTime(j,iDetsIn)'-min(TrigTime(j,iDetsIn));
        if CORREL
            correlTime = TrigCor(j,iDetsIn)';
            correlCoef = CoefCor(j,iDetsIn)';
        else
            correlTime = zeros(nDetsIn,1);
            correlCoef = zeros(nDetsIn,1);
        end
        Sat = (MaxRaw(j,iDetsIn)==255 | MinRaw(j,iDetsIn)==0)';
        Amp = AmpMax(j,iDetsIn)';
        txtTable(1:nDetsIn,:) = [uTime AntId EvtId CoincId trigTime correlTime correlCoef Sat Amp zeros(nDetsIn,1)];
        txtTable = sortrows(txtTable,1);
        
        for k = 1:size( txtTable, 1 )
            fprintf( fid, '%20f ', txtTable( k, 1 ) );   % Event time after correction (in sample counts)
            fprintf( fid, '%6d ',  txtTable( k, 2) );    % Antenna ID 
            fprintf( fid, '%6d ',  txtTable( k, 3) );    % Event number
            fprintf( fid, '%6d ',  txtTable( k, 4) );    % Coinc number
            fprintf( fid, '%8.2f ', txtTable( k, 5 ) );  % Delay with 1rst antenna (standard method)
            fprintf( fid, '%8.2f ', txtTable( k, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
            fprintf( fid, '%6.2f ', txtTable( k, 7 ) );  % Average correlation for this antenna (intercorrelation method)
            fprintf( fid, '%6d ', txtTable( k, 8 ) );  %  Saturation flag
            fprintf( fid, '%6.2f ', txtTable( k, 9 ) );  % Amplitude
            fprintf( fid, '%6.3f ', 0 ); % Gain OBSOLETE
            fprintf( fid, '%6.3f ', 0);  % Amplitude method 1 (bline)
            fprintf( fid, '%6.3f ', 0);  % Amplitude method 2 (PSD)
            fprintf( fid, '\n' );
        end
        clear txtTable
        
        decim = floor(length(ind_val)/10);
        if floor(j/decim)==j/decim
            disp(sprintf('Done at %2.0f percent',j/decim*10));
        end
    end
end

if FlagTable==0
    if Type==0
        disp(sprintf('No antenna coincs in R%d!',nrun));
    elseif Type==1
        disp(sprintf('No scintillator coincs in R%d!',nrun));
    elseif Type==2
        disp(sprintf('No hybrid coincs in R%d...',nrun));
    end
    return
end


%% Output to file
txtTable(1,:)=zeros(1,10);  % !!!! Hard cod� !!!!!
txtTable(1,2) = Detectors(1);
for k = 1:size( txtTable, 1 )
	fprintf( fid, '%20f ', txtTable( k, 1 ) );   % Event time after correction (in sample counts)
	fprintf( fid, '%6d ',  txtTable( k, 2) );    % Antenna ID 
	fprintf( fid, '%6d ',  txtTable( k, 3) );    % Event number
	fprintf( fid, '%6d ',  txtTable( k, 4) );    % Coinc number
	fprintf( fid, '%8.2f ', txtTable( k, 5 ) );  % Delay with 1rst antenna (standard method)
	fprintf( fid, '%8.2f ', txtTable( k, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
	fprintf( fid, '%6.2f ', txtTable( k, 7 ) );  % Average correlation for this antenna (intercorrelation method)
	fprintf( fid, '%6d ', txtTable( k, 8 ) );  %  Saturation flag
	fprintf( fid, '%6.2f ', txtTable( k, 9 ) );  % Amplitude
	fprintf( fid, '%6.3f ', txtTable( k, 10 ) );  % Gain
	fprintf( fid, '%6.3f ', 0);  % Gain from PSD (method 2)
	fprintf( fid, '%6.3f ', 0);  % Gain from PSD (method 1)
	fprintf( fid, '\n' );
end
clear txtTable

disp(sprintf('%d coincs written to R%d_coinctable.txt file.',n,nrun));
fclose( fid );

