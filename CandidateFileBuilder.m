function [Struct]=CandidateFileBuilder(Struct)
% Builds table for C++ shower reconstruction 
% 18/05/2012 OMH

fclose all;

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
%% Change to first 4 numbers only (in case of meta run)
nrun_temp=num2str(nrun);
nrunpsd=str2num(nrun_temp(1:4));

SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
CoincStruct = Struct.Coinc;
CoincId = CoincStruct.IdCoinc;
tag = CoincStruct.Det.Tag;
sig = CoincStruct.Det.Sigma;
rawamp = CoincStruct.Det.AmpMax;
Status = CoincStruct.Det.Status;
if ~isfield(Struct.Coinc,'DelayCorrRecons')
    thetas = CoincStruct.SphRecons.Theta;
    flags = CoincStruct.SphRecons.Flag;
    thetap = CoincStruct.PlanRecons.Radio.Theta;
    flagp = CoincStruct.PlanRecons.Radio.Flag;
    MinDist = CoincStruct.SphRecons.minDistSource;
    slopep = CoincStruct.PlanRecons.Radio.SlopeDelay;
    chi2p = CoincStruct.PlanRecons.Radio.Chi2Delay;
    slopes = CoincStruct.SphRecons.SlopeDelay;
    chi2s = CoincStruct.SphRecons.Chi2Delay;
else
    thetas = CoincStruct.DelayCorrRecons.SphRecons.Theta;
    flags = CoincStruct.DelayCorrRecons.SphRecons.Flag;
    thetap = CoincStruct.DelayCorrRecons.PlanRecons.Radio.Theta;
    flagp = CoincStruct.DelayCorrRecons.PlanRecons.Radio.Flag;
    MinDist = CoincStruct.DelayCorrRecons.SphRecons.minDistSource;
    slopep = CoincStruct.DelayCorrRecons.PlanRecons.Radio.SlopeDelay;
    chi2p = CoincStruct.DelayCorrRecons.PlanRecons.Radio.Chi2Delay;
    slopes = CoincStruct.DelayCorrRecons.SphRecons.SlopeDelay;
    chi2s = CoincStruct.DelayCorrRecons.SphRecons.Chi2Delay;  
end;

%
doShowerRecons = zeros(1,length(CoincId));
ShowerSignals = cell(length(CoincId),length(Detectors));

%% Load calib DST
dstname = [CAL_PATH sprintf(calib_filename,nrunpsd,1)];
if fopen(dstname)>0
    disp(sprintf('Calib DST %s found.',dstname));  
    calibdst = load(dstname);
    antcal = [calibdst.Struct.Setup.Antennas];
    nocal = calibdst.Struct.Calib.NoCalib;
    calant = calibdst.Struct.Setup.Antennas;
    if size(setdiff(calant,nocal),2)==0
         disp 'No antennas calibrated in this DST!'
         calibin = 0;
    else
         calibin = 1;     
         gaincal = [calibdst.Struct.Calib.GainLin];
         [a ia iantcal] = intersect(Detectors,antcal);
    end
else
    disp(sprintf('Calib DST %s not found.',dstname));  
    calibin = 0;
end

%% Dump coordinates   
DumpCoordTxt(Struct);


%% Compute Amplitudes & load signals to structure
CalAmp1 = zeros(length(CoincStruct.IdCoinc),length(Detectors));
CalAmp2 = zeros(length(CoincStruct.IdCoinc),length(Detectors));
for i = 1:length(CoincId)  % loop on all valid coincs
    iDetsIn =  find(tag(i,:)==1);  % Triggering detecors
    iAntsIn =  find(tag(i,:)==1 & DetectorType==0);  % Triggering antennas
    Antennas = Detectors(iAntsIn);
    nAntsIn = size(iAntsIn,2);
    trig = zeros(1,length(Detectors));
    for j = 1:length(Detectors)
      st = dec2bin(Status(i,j));
      stv = int16(sscanf(st,'%1d'))'; %LSB last
      stv=fliplr(stv); %LSB first
      trig(j) = stv(1);
    end
    
    SelectionTag = nAntsIn>nAntsCut & MinDist(i)>MinDistCut & thetap(i)<=ThetaCut & thetas(i)<=ThetaCut & flags(i)==1 & flagp(i)==1 & chi2s(i)<Chi2Cut & chi2p(i)<Chi2Cut & abs(1-slopes(i))<SlopeCut & abs(1-slopep(i))<SlopeCut ;

    % Here we select coincs with valid criteria for shower recons
    %display(sprintf('Nbant: %d, Nbsci: %d, ThetaP: %d, FlagP: %d, ThetaS: %d, FlagS: %d',nAntsIn,nbsci,thetap(j),flagp(j),thetas(j),flags(j)))
    if (SelectedData | (CandidateSelection==1 & SelectionTag))
        doShowerRecons(i) = 1;   
        
        % Method 1: use Sigma as calib coef
        CalAmp1(i,iAntsIn) = rawamp(i,iAntsIn)./sig(i,iAntsIn);
        
        % Method 2: use Calib DST
        if calibin 
          aout = intersect(Antennas,nocal);
          if size(aout,2)==0  % DST calib available for all antennas in candidate
            for j = 1:length(iAntsIn)
                ant = Antennas(j);
                calantind = find(antcal==ant);
                CalAmp2(i,iAntsIn(j)) = rawamp(i,iAntsIn(j))./gaincal(calantind)*1e6; % in muV
            end
            [m indm] = max(CalAmp2(i,:));
            CalAmp1(i,:) = CalAmp1(i,:)*CalAmp2(i,indm)./CalAmp1(i,indm);
          else
            disp(sprintf('Could not find calib DST for Antenna %d at location %s.',aout(1),dstname));  
            CalAmp2(i,iAntsIn) = -1;
          end  
        else
            CalAmp2(i,iAntsIn) = -1;
        end
        %[Antennas' CalAmp1(i,iAntsIn)' CalAmp2(i,iAntsIn)']
        
        % Now load signals
        trigIn = find(trig == 1);
        Dets = Detectors(trigIn);
        Events = CoincStruct.Det.Evt(i,trigIn);
        for j=1:length(Dets)
          %disp(sprintf('Event %d Detector %d',Events(j),Dets(j)))
          % Get data
          fd = OpenFileData( nrun, Dets( j ) );
          if fd<0
            disp(sprintf('Coinc %d: could not find data for detector %d.',CoincId(i), Dets(j)))
            pause
            fclose all;
            continue
          end
          fseek( fd, ibuff*(Events(j)-1),'bof');
          DataEvt = double( fread( fd, ibuff, 'uint8' ) );
%           figure; plot(DataEvt)
%           pause
          ShowerSignals{ i, trigIn(j) } = DataEvt*SCALE; % Now in Volts    
          fclose(fd);
        end
    end
end

if ~isfield(Struct.Coinc,'DelayCorrRecons')
    Struct.Coinc.IsShower = doShowerRecons;
else
    Struct.Coinc.DelayCorrRecons.IsShower = doShowerRecons;
end;
Struct.Coinc.Det.CalibratedAmp1 = CalAmp1;
Struct.Coinc.Det.CalibratedAmp2 = CalAmp2;

if RecordSignals==1
    if ~isfield(Struct.Coinc,'DelayCorrRecons')
        Struct.Coinc.ShowerRecons.ShowerSignals = ShowerSignals;
    else
        Struct.Coinc.DelayCorrRecons.ShowerRecons.ShowerSignals = ShowerSignals;
    end
else
    if ~isfield(Struct.Coinc,'DelayCorrRecons')
        Struct.Coinc.ShowerRecons.ShowerSignals = cell(length(CoincId),length(Detectors));
    else
        Struct.Coinc.DelayCorrRecons.ShowerRecons.ShowerSignals = cell(length(CoincId),length(Detectors));
    end
end;

%% Load & complete coinctable
filename = [TEXT_PATH sprintf( 'R%d_coinctable.txt', nrun)];
if fopen( filename )>0
    fclose all;
    coincTable = load( filename );   
    if size(coincTable,1) == 0
       fclose all
       disp(sprintf('File %s is empty.',filename));
       return
    end
else
    fclose all;
    disp(sprintf('File %s does not exist.',filename));
    Struct.Coinc.ShowerRecons.ShowerSignals = ShowerSignals;
    return
end  
IdCoincs = coincTable(:,4);
coincTable(:,11)=0;  % Add column for bline gain
coincTable(:,12)=0;  % Add column for PSD gain


% Write results to coinctable & structure
ShowerId = CoincId(find(doShowerRecons==1));


if sum(doShowerRecons)==0
    disp(sprintf('No shower candidates in R%d.',nrun));
else
    disp(sprintf('%d shower candidates in R%d.',sum(doShowerRecons),nrun));    
    for i = 1 : length(ShowerId)
        
        % Now reorder calibrated amplitudes following increasing time
        ind = find(IdCoincs == ShowerId(i));  % Index for txt table
        ind_struct = find(CoincId == ShowerId(i)); % Index for dst
        
        order = zeros(1,length(ind));
        for k=1:length(ind)
            order(k)=find(Detectors==coincTable(ind(k),2));
        end
        coincTable(ind,11) = CalAmp1(ind_struct,order);
        coincTable(ind,12) = CalAmp2(ind_struct,order);  
    end
end



%[coincTable( :, 4), coincTable(:,2), coincTable(:,11), coincTable(:,12)] 


%% Update coinctable file
fclose all;
filename = [TEXT_PATH sprintf( 'R%d_coinctable.txt', nrun)];
fid = fopen( filename, 'w' );
for k = 1:size( coincTable, 1 )
  fprintf( fid, '%20f ', coincTable( k, 1 ) );   % Event time after correction (in sample counts)
  fprintf( fid, '%6d ',  coincTable( k, 2) );    % Antenna ID 
  fprintf( fid, '%6d ',  coincTable( k, 3) );    % Event number
  fprintf( fid, '%6d ',  coincTable( k, 4) );    % Coinc number
  fprintf( fid, '%8.2f ', coincTable( k, 5 ) );  % Delay with 1rst antenna (standard method)
  fprintf( fid, '%8.2f ', coincTable( k, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
  fprintf( fid, '%6.2f ', coincTable( k, 7 ) );  % Average correlation for this antenna (intercorrelation method)
  fprintf( fid, '%6d ', coincTable( k, 8 ) );  %  Saturation flag
  fprintf( fid, '%6.3f ', coincTable( k, 9 ) );  % Amplitude
  fprintf( fid, '%6.3f ', coincTable( k, 10 ) );  % Gain
  fprintf( fid, '%8.2f ', coincTable( k, 11 ) );  % Gain from baseline 
  fprintf( fid, '%8.2f ', coincTable( k, 12 ) );  % Gain from PSD 
  fprintf( fid, '\n' );
end
fclose( fid );
clear coincTable
