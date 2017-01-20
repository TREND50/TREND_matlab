function [Struct]=CalibStructureBuilder(nrun)
% Build calibration structure
% based on PSD files only.
% OMH 28/02/2012

nrun = sscanf(nrun,'%d')

% load global variables
SharedGlobals;
BWmin = 55e6;
BWmax = 95e6;
thresh = 0.4; %gain varation/min
Npsd = 20; % Number of DST runs to scan
difthr = 0.1;  % 10%

%% Setup
setup.Run = nrun;
setup.Antennas = [101:138 140 148:158]; % Hardcoded
nant = length(setup.Antennas);
setup.LoadRun = zeros(1,nant);
setup.RefRun = zeros(1,nant);

if nrun<2900
  setup.LoadRun(1:nant) = 2667;  % Run performed with load to determine the gain of the acquisition chain
  setup.RefRun(1:nant) = 2676;  % Run performed with antenna just after to serve a s a reference for the gain evolution 
  slicestart = 2600;
  slicestop = 2900;
elseif nrun>=2900 & nrun<3087
  setup.LoadRun(1:nant) = 2987;  % Run performed with load to determine the gain of the acquisition chain
  setup.RefRun(1:nant) = 2988;  % Run performed with antenna just after to serve a s a reference for the gain evolution 
  setup.RefRun(setup.Antennas>=148) = 3004;
  setup.LoadRun(setup.Antennas==157) = 3129;
  setup.RefRun(setup.Antennas==157) = 3130;
  setup.LoadRun(setup.Antennas==158) = 3129;
  setup.RefRun(setup.Antennas==158) = 3130;
  slicestart = 2900;
  slicestop = 3087;
elseif nrun>=3087 & nrun<3775
%   Calib performed in October 2011 (3129&3130) are not reliable because
%   antennas were cut afterwards.
%   setup.LoadRun(1:nant) = 3129;  % Run performed with load to determine the gain of the acquisition chain
%   setup.RefRun(1:nant) = 3130;  % Run performed with antenna just after to serve a s a reference for the gain evolution 
%   setup.LoadRun(setup.Antennas==106) = 2987;
%   setup.RefRun(setup.Antennas==106) = 2988;
%   setup.LoadRun(setup.Antennas==109) = 2987;
%   setup.RefRun(setup.Antennas==109) = 2988;
%   setup.RefRun(setup.Antennas==112) = 3155;
%   slicestart = 3087;
%   slicestop = 3200;
  setup.LoadRun(1:nant) = 3759;  % Run performed with load to determine the gain of the acquisition chain
  setup.RefRun(1:nant) = 3761;  % Run performed with antenna just after to serve a s a reference for the gain evolution 
  slicestart = 3087;
  slicestop = 3775;
elseif nrun>=3775
  setup.LoadRun(1:nant) = 3793;  % Run performed with load to determine the gain of the acquisition chain
  setup.RefRun(1:nant) = 3794;  % Run performed with antenna just after to serve a s a reference for the gain evolution 
  slicestart = 3775;
  slicestop = 4442;
end
%
setup.LoadRunDate = zeros(1,nant);
setup.RefRunLST = zeros(1,nant);
setup.RefRunDate = zeros(1,nant);
setup.PresRun = zeros(1,nant);
setup.PresRunDate = zeros(1,nant);
setup.PresRunLST = zeros(1,nant);
setup.PresRunInd = zeros(1,nant);
setup.TimeStep = zeros(1,nant);
%
calib.G0dB = zeros(1,nant);
calib.meanPSDLoad = zeros(1,nant);
calib.meanPSDRef = zeros(1,nant);
calib.meanPSDPres = zeros(1,nant);
calib.GainLin = zeros(1,nant);
calib.Flag = zeros(1,nant);
calib.OptCalibCoefs = zeros(4,nant);

%% Load optical system calib curves
% Old systems 
if nrun<3300
    filename = [CALIB_PATH 'CalibOptical_old.txt'];
    if fopen(filename)
        optcalib = load(filename);
        antid = optcalib(:,1);
        parCal = optcalib(:,2:5);
    end
else
    antid = [];
    parCal = [];
end

%% % Load current psd & compute status
d = ReadPSD( setup.Run, setup.Antennas );
if size(d,1)==0
    disp(sprintf('Could not find current PSD run (R%d). Aborting.',setup.Run ))
    Struct.Setup = setup;
    Struct.Calib = calib;    
    return
end
setup.RunDate = max([d(:).tdeb]);
if size(setup.RunDate,1)==0
    setup.RunDate = 0;
end

Nants = size(d,2); % Number of antennas with gains computed
msiz = 0;
for i=1:Nants
  msiz = max(msiz,length(d(i).t));
end
calib.Status = -ones(nant,msiz);
calib.meanPSDEvo = zeros(nant,msiz);
nscan = zeros(1,Nants);

for i=1:Nants % loop on antennas
  if length(d(i).t) == 0
      disp(sprintf('Error! Size of time vector is %d for antenna %d in PSD%d.',length(d(i).t),d(i).id,setup.Run))
      continue;
  end
  ind = find(setup.Antennas==d(i).id);
  calib.Time(ind,1:length(d(i).t)) = d(i).t;  
  if size(ind,2) == 0
      disp(sprintf('Error! Antenna %d in PSD file.',d(i).id))
      return;
  end
  % Compute status
  psd_t = d(i).t;
  stept = mean(diff(psd_t))/60;  % Minutes
  steph = stept/60;  % Hours
  setup.TimeStep(i) = round(stept);
  bit0 = zeros(1,length(psd_t));
  bit1 = zeros(1,length(psd_t));
  bit2 = zeros(1,length(psd_t));
  bit3 = zeros(1,length(psd_t));
  bit4 = zeros(1,length(psd_t));
  
  for k = 1:length(psd_t)  % loop on time
    CurrentPSD = d( i ).psd( :, k ); % Current PSD
    ref = mean( CurrentPSD(d(i).f > 10e6 & d(i).f < 20e6));
    mg(k) = mean( CurrentPSD(d(i).f > BWmin & d(i).f < BWmax)); % Average PSD for the reference run
    calib.meanPSDEvo(ind,k) = mg(k); 
    flat(k) = std( CurrentPSD(d(i).f > BWmin & d(i).f < BWmax)-mg(k) );
    if mg(k)-ref<4
        bit0(k) = 1; % Antenna is dead (gain extremely low)
    end
    if flat(k)>2
        bit1(k) = 1; % Bad shape for PSD
    end
    if k>3 & abs(mg(k-2)-mg(k-1))>thresh*stept & abs(mg(k-3)-mg(k-2))<thresh*stept & abs(mg(k-1)-mg(k))<thresh*stept & flat(k-2)<2 & flat(k-1)<2            
        if mg(k-1)>mg(k-2)
            bit2(k-2) = 1; %Switch to high state
        else
            bit3(k-2) = 1;    %Switch to low state
        end
        calib.Flag(ind) = 3;   % Jump
    end
    nscan(k) = ceil(0.5/steph);
    if k>2*nscan(k)+1 
        selmg = mg(k-(2*nscan(k)+1):k);
        out = find(abs(mg(k-nscan(k)-1)-selmg)>1);
        if length(out)>1
            bit4(k-nscan(k)-1) = 1;  % More than 1 point with variation of PSD>1dB within +-1h around measurement
        end
    end
    if DISPLAY
      setup.Antennas(i)
      figure(1)    
      plot(d(i).f ,CurrentPSD,'k')
      grid on
      hold on
      line([BWmin BWmax],[mg(k) mg(k)])
    end
  end
  calib.Status(ind,1:length(psd_t)) = bit4.*2^4 + bit3.*2^3 + bit2.*2^2 + bit1.*2^1 + bit0.*2^0;
end

%% Get calib runs (load+reference run)

%% Compute calib coefficient
% Preamp thermal noise
kB   = 1.38e-23; % Boltzmann
T    = 273 + 15; % Temperature at measurment time (�K)
R    = 75; % Load impedance
NF   = 0.4;  % LNA noise factor
psd_in = 10*log10( kB*T*R ) + NF;  % PSD at LNA input
%

for i=1:nant
  % Load load run data
  d = ReadPSD( setup.LoadRun(i), setup.Antennas(i) );
  if size(d,1)==0
    disp(sprintf('Could not find PSD reference run %d for antenna %d. Skip.',setup.RefRun(i),setup.Antennas(i) ))
    continue;
  end
  ind = find(setup.Antennas==d.id);
  setup.LoadRunDate(ind) = d.tdeb; 
  lastpsd = d.psd(:,end);  % Take last PSD measurment
  if setup.Antennas(i)  == 125 & setup.LoadRun(i)==2987  % Patch for antenna1265 & Load Run = 2987 (gain drops during run)
      lastpsd = d.psd(:,1);  % Take 1st PSD measurment
  end
  calib.meanPSDLoad(ind) = mean(lastpsd( d.f>BWmin & d.f<BWmax));
  calib.G0dB(ind) = calib.meanPSDLoad(ind) - psd_in;
end

%

for i=1:nant
  % Load reference antenna run data
  d = ReadPSD( setup.RefRun(i), setup.Antennas(i) );
  if size(d,1)==0
    disp(sprintf('Could not find PSD reference run %d for antenna %d. Skip.',setup.RefRun(i),setup.Antennas(i) ))
    continue;
  end
  ind = find(setup.Antennas==d.id);
  setup.RefRunDate(ind) = d.tdeb;
  [year, month, day, hour, minute, second] = UnixSecs2Date(setup.RefRunDate(ind));
  setup.RefRunLST(ind) = utdate2lst([year, month, day, hour, minute, second]);
  firstpsd   = d.psd(:,1)';  % Average PSD for the reference run
  calib.meanPSDRef(ind) = mean( firstpsd( d.f > BWmin & d.f < BWmax)); % Average PSD for the reference run
end


%% Get present PSD level and compute gain
getDailyPSD = zeros(nant,1);
%LST time of ref run... DailyPSD has to be performed at a time with similar level (+- 5%)
% Alternative would be to get ALL times with simiilar Galactic noise... 
GalTime = [0:0.5:24];
GalNoise = [ 1.4556 1.4482 1.4373 1.4293 1.4256 1.4238 1.4207 1.4164 1.4084 1.3951 1.3757 1.3485 1.3134 1.2702 1.2224 1.1721 1.1214...
             1.0770 1.0395 1.0145 1.0000 1.0008 1.0146 1.0407 1.0766 1.1225 1.1767 1.2436 1.3254 1.4263 1.5481 1.6847 1.8266 1.9614...
             2.0784 2.1633 2.2072 2.2086 2.1736 2.1038 2.0140 1.9109 1.8075 1.7143 1.6332 1.5701 1.5225 1.4886 1.4556];  % Theoeretical Galactic noise (normalised)
% figure;
% plot(GalTime,GalNoise,'k')
% hold on
% plot(GalTimeMod,GalNoiseMod,'r')
%GalNoiseRefRun = GalNoise(round(lstRefRun*2));  % Galactic Noise level at time of RefRun.
%GalNoiseRefInd = find(abs(GalNoise-GalNoiseRefRun)/GalNoiseRefRun<0.05)/2 % Time window to seek for present PSD measurment (ie same -5% level- Galactic conditions) 

% Loop on runs while a calib point could not be be computed for all of them
i = 0;
ncount = 0;
while sum(getDailyPSD)<nant & ncount<=Npsd
    runpsd = nrun+(-1)^mod(i,2)*round(i/2); % look once later, once sooner
    if runpsd>=slicestop | runpsd<slicestart
        i = i+1;
        continue % Scan PSD run in the same time slice as nrun (suffisant condition for same setup).
    end
    disp(sprintf('Looking for PSD run %d...',runpsd));
    d = ReadPSD(runpsd,setup.Antennas);
    if size(d,1)>0  % PSD data was found in this run 
       Nmes = size(d(1).t,2);       
       Nants = size(d,2); % Number of antennas with PSD in this run
       disp(sprintf('Data found for PSD run %d: %d antennas and %d measurments',runpsd,Nants,Nmes));
       if d(1).t(end)<1  % Run is less than 1s.... Used to be 3600s
           disp(sprintf('PSD run is only %3.1f h. Too short to compute gain. Skip',d(1).t(end)/3600))
           i = i+1;
           continue
       end
       % Look for Calib dst       
       dstfile = [CALIB_PATH sprintf(calib_filename,runpsd,1)];
       %if 0
       if runpsd<nrun & fopen(dstfile)>0
           disp(sprintf('Calib dst found for R%d.',runpsd));
           calibdst = load(dstfile);
           status = calibdst.Struct.Calib.Status;
           steph = setup.TimeStep/60;
           for j=1:Nants % loop on antennas
              ind = find(setup.Antennas==d(j).id); 
              nscan(j) = ceil(0.5/steph(ind));  
           end
       else  % No calib dst for this run. Status has to be computed...
          status = -ones(Nants,Nmes);
          for j=1:Nants % loop on antennas
              ind = find(setup.Antennas==d(j).id);  
              if size(ind,2) == 0
                  disp(sprintf('Error! Antenna %d in PSD file.',d(j).id))
                  return;
              end
              psd_t = d(j).t;
              stept = mean(diff(psd_t))/60; % min
              steph = stept/60; % hour
              bit0 = zeros(1,length(psd_t));
              bit1 = zeros(1,length(psd_t));
              bit2 = zeros(1,length(psd_t));
              bit3 = zeros(1,length(psd_t));
              bit4 = zeros(1,length(psd_t));
              
              for k = 1:length(psd_t)  % loop on time
                CurrentPSD = d( j ).psd( :, k ); % Current PSD
                ref = mean( CurrentPSD(d(j).f > 10e6 & d(j).f < 20e6));
                mg(k) = mean( CurrentPSD(d(j).f > BWmin & d(j).f < BWmax)); % Average PSD for the reference run
                flat(k) = std( CurrentPSD(d(j).f > BWmin & d(j).f < BWmax)-mg(k) );
                % From most to least serious (following elseif not considered as soon as one is true)
                if mg(k)-ref<4
                    bit0(k) = 1; % Antenna is dead (gain extremely low)
                end
                if flat(k)>2
                    bit1(k) = 1; % Bad shape for PSD
                end
                if k>3 & abs(mg(k-2)-mg(k-1))>thresh*stept & abs(mg(k-3)-mg(k-2))<thresh*stept & abs(mg(k-1)-mg(k))<thresh*stept & flat(k-2)<2 & flat(k-1)<2
                    if mg(k-1)>mg(k-2)
                      bit2(k-2) = 1; %Switch to high state
                    else
                      bit3(k-2) = 1;    %Switch to low state
                    end
                end
                nscan(j) = ceil(0.5/steph);
                if k>2*nscan(j)+1 
                    selmg = mg(k-(2*nscan(j)+1):k);
                    out = find(abs(mg(k-nscan(j)-1)-selmg)>1);
                    if length(out)>1
                        bit4(k-nscan(j)-1) = 1;  % Variation of PSD>1dB within +-1h around measurement
                    end
                end
              end
              status(ind,1:length(psd_t)) = bit4.*2^4 + bit3.*2^3 + bit2.*2^2 + bit1.*2^1 + bit0.*2^0;
          end
       end
       for k=1:Nants   % loop on antennas
         ind = find(setup.Antennas==d(k).id);
         if getDailyPSD(ind)==1  % Daily PSD already computed for this antenna. Skip it  
             continue
         end
         disp(sprintf('*** Computing calib coef for antenna %d ***', setup.Antennas(ind)));
                
         sizt = length(d(k).t);
         dstlsttime = zeros(1,sizt);
         good = zeros(1,sizt);
         good( status(ind,1:sizt) == 0) = 1;
         good(1:min(sizt,nscan(k))) = 0; % Skip nscan 1st measurments because jump could not be assessed 
         good(max(1,sizt-nscan(k)):sizt) = 0;  % Skip nscan last measurments because jump could not be assessed 
         [year, month, day, hour, minute, second] = UnixSecs2Date(max([d(:).tdeb])+d(k).t);
         selgood = find(good==1);
         if size(selgood,2) == 0
             disp(sprintf('1: No valid measurment in PSD%d for antenna %d',runpsd,setup.Antennas(ind)))  
             continue
         end
         for j=1:length(selgood)   
              % Get time of measurement closest to that of Ref run...
              indj = selgood(j);
              dstlsttime(indj) = utdate2lst([year(indj), month(indj), day(indj), hour(indj), minute(indj), second(indj)]);
         end
         %
         dstlsttime(good==0)=-1e12;  % Get rid of these...
         [a indt] = min(abs(dstlsttime-setup.RefRunLST(ind)));
         dstlsttime(good==0)=0;  % Put them back  
         
         % Compare the Gal envir at measurement time and Ref run time.
         GalNoise_at_RefTime = GalNoise(max(1,round(2*setup.RefRunLST(ind))));
         GalNoise_at_PresentTime = GalNoise(max(1,round(2*dstlsttime(indt ))));
         dif = abs(GalNoise_at_PresentTime-GalNoise_at_RefTime)/GalNoise_at_RefTime;
         if size(dif, 2)>0
           disp(sprintf('Ref measurement: %3.3f at %3.1f h LST',GalNoise_at_RefTime, setup.RefRunLST(ind)))
           disp(sprintf('Best match: %3.3f at %3.1fh LST. Dif = %3.1f pc ',GalNoise_at_PresentTime,dstlsttime(indt),dif*100))
         else
           disp(sprintf('2: No valid measurment in PSD%d for antenna %d',runpsd,setup.Antennas(ind)))  
         end
         %
         if size(dif, 2)>0 & dif<difthr % Gal envir similar at the "difthr" level -> compute gain
            % Check current PSD quality           
            setup.PresRun(ind) = runpsd;
            setup.PresRunDate(ind) = max([d(:).tdeb]);
            setup.PresRunLST(ind) = dstlsttime(indt);
            setup.PresRunInd(ind) = indt;
            CurrentPSD = d(k).psd( :, indt ); % Current PSD
            calib.meanPSDPres(ind) = mean( CurrentPSD(d(k).f > BWmin & d(k).f < BWmax)); % Average PSD for the reference run at the reference time
            GaindB = calib.G0dB(ind) + calib.meanPSDPres(ind) - calib.meanPSDRef(ind);
            calib.GainLin(ind) = 10.^(GaindB./20);
            getDailyPSD(ind)=1; 
            % Load optical system calibration curve
            indcal = find(antid==setup.Antennas(ind));
            if size(indcal,1)>0
                calib.OptCalibCoefs(:,ind) = parCal(indcal,:);
            end
            disp(sprintf('PSD%d valid for calibration for antenna %d in R%d : GainLin = %3.2e',setup.PresRun(ind),setup.Antennas(ind), nrun, calib.GainLin(ind)))
         end
       end
       ncount = ncount+1;
    end
    i = i+1;   
    fclose all;
end
if sum(getDailyPSD)<nant
    disp('Warning: no gain could be computed for the following antennas:')
    setup.Antennas(getDailyPSD==0)
    calib.NoCalib = setup.Antennas(getDailyPSD==0);
    
end

%% Saving DST
Struct.Setup = setup;
Struct.Calib = calib;
filename = [DST_PATH sprintf(calib_filename,nrun,1)];
save(filename,'Struct');

       


    
    
