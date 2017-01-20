function [Struct]=StatusBuilder(Struct)
% Complete status info
% 13/06/11 OMH
% Added info from calib dst
% 09/05/12 OMH

%% Initialisation
RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;        
Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
Evt=[RunSetup.Det.Evt];
CoincStruct = Struct.Coinc;
minraw = CoincStruct.Det.MinRaw;
maxraw = CoincStruct.Det.MaxRaw;
status = CoincStruct.Det.Status;
IdCoinc = CoincStruct.IdCoinc;
unixs = CoincStruct.Det.UnixTime;
trigRate = [RunSetup.InfosRun.TrigRate];
trigTime = [RunSetup.InfosRun.TrigTime];
NbCoinc = length(IdCoinc);
calibin = 0;

%% Change to first 4 numbers only (in case of meta run)
nrun_temp=num2str(nrun);
nrunpsd=str2num(nrun_temp(1:4));
%% Load calib DST
dstname = [CAL_PATH sprintf(calib_filename,nrunpsd,1)];
if fopen(dstname)>0
     disp(sprintf('Calib DST %s found.',dstname));  
     calibdst = load(dstname);
     nocalant = calibdst.Struct.Calib.NoCalib;
     calant = calibdst.Struct.Setup.Antennas;
     if size(setdiff(calant,nocalant),2)==0
         disp 'No antennas calibrated in this DST!'
         calibin = 0;
     else
         calibin = 1;     
         utime_cal = calibdst.Struct.Calib.Time+calibdst.Struct.Setup.RunDate; % Time of PSD measurments ins seconds
         status_cal = calibdst.Struct.Calib.Status;
     end
else
     disp(sprintf('Calib DST %s not found.',dstname));  
end

    %% Loop on coincs
 for i = 1:NbCoinc
    utime = max(unixs(i,:));
%    disp(sprintf('Coinc %d at %d secs',IdCoinc(i),utime));
    ind = find(trigTime<utime,1,'last'); % Get time period of trigger rate measurement coincident with coincidence
    for j = 1:length(Detectors)        
        %% High trigger rate flag
        trate = trigRate(ind,j); 
%        if ind<length(trigTime)
%              disp(sprintf('Average trigger rate on antenna %d in the time period %d-%d s = %3.1f Hz',Detectors(j),trigTime(ind),trigTime(ind+1),trate))
%        else
%              disp(sprintf('Average trigger rate on antenna %d in the time period > %d s = %3.1f Hz',Detectors(j),trigTime(ind),trate))
%        end
        if trate>TrigRateLimit
            status(i,j) = status(i,j) + 4;  % Bit2 set at 1
        end
        
        %% Saturation flag
        if (minraw(i,j) == 0 | maxraw(i,j) == 2^NBITS)
            status(i,j) = status(i,j) + 8; % Bit3 set at 1
        end
        
        %% PSD status (calib DST)
        if calibin
          indcal = find(calibdst.Struct.Setup.Antennas==Detectors(j));  
          if size(indcal,2)>0
            if calibdst.Struct.Calib.Flag(indcal)>0
              status(i,j) = status(i,j)+256; % Bit8 set at 1 if Jump flag
            end
            [a indt] = min(abs(utime_cal(indcal,:)-utime));
            if (status_cal(indcal,indt)>=0) % DST data is available for this anenna at this time
                status(i,j) = status(i,j)+16; % Bit4 set at 1 if Calib DST data available
                st = dec2bin(status_cal(indcal,indt));
                %Detectors(j)
                %status_cal(indcal,indt)
                stv = int16(sscanf(st,'%1d'))'; %LSB last
                if size(stv,2) == 5
                  stv(2:3) = []; % Remove jump tags
                elseif size(stv,2) == 4
                  stv(1:2) =  []; % Remove jump tag
                elseif size(stv,2) == 3
                  stv(1) =  []; % Remove jump tag
                end
                stv = [stv 0 0 0 0 0];  % Add the 5 fields corresponding to the event itself
                %stv=fliplr(stv);
                pow = numel(stv)-1:-1:0;
                exp = int16(2.^pow);
                stat_psd = sum(stv.*exp); % Convert bit pattern to decimal value
                status(i,j) = status(i,j)+stat_psd;
            end
          end
        end
    end
end
Struct.Coinc.Det.Status = status;
