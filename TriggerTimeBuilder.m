function [Struct]=TriggerTimeBuilder(Struct)
% Compute trigger times
% 06/12/10, OMH
% Last modification : 


RunSetup = Struct.Setup;
nrun = RunSetup.Run;
if nrun < 1992
    disp 'Soft not valid for this run (pulse has to be centered)'.
    return
end

SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
Evt=[RunSetup.Det.Evt];
NbCoinc=RunSetup.TotalCoinc;

CoincStruct = Struct.Coinc;
tag = CoincStruct.Det.Tag;

AmpMax = zeros(NbCoinc,length(Detectors));
TrigTime = zeros(NbCoinc,length(Detectors));
Gain = zeros(NbCoinc,length(Detectors));

%% Loop on detectors
for i=1:length(Detectors)
    if exist('id') & Detectors(i)~=id
        continue
    end
    if DetectorType(i) == 0 % antennas
        ib = ibuff;
    else
        ib = ibuffs;
    end
    ts = [ 1: ib  ]/FSAMPLING;
    fover = 10;
    tsover = [1 :1/fover: ib]/FSAMPLING;
    %tmu = ts*1e6;
    
    % Open data file
    fd=OpenFileData(nrun,Detectors(i));
    if fd<0
        disp(sprintf('Data file not found for detector %d in run %d!',Detectors(i),nrun));
        pause
        continue
    end
    if ismember(nrun,runs2010) % Old soft
        ft=OpenFileTime(nrun,Detectors(i));
    end
    ind_tag=find(tag(:,i)==1);
    
    %% Loop on events
    for j=1:length(ind_tag)
        deltat = 10;
        evt = CoincStruct.Det.Evt(ind_tag(j),i);
        t = CoincStruct.Det.Time(ind_tag(j),i);        
        in = find(tag(ind_tag(j),:)==1);  % Events with good quality (==2) or triggering only (==1)
        
        if exist('event') & event~=evt
            continue
        end
        fseek(fd,ib*(evt-1),'bof');
        DataEvt=double(fread(fd,ib,'uint8'));
                
        if length(DataEvt)~=ib
            disp(sprintf('Error! Data size for event %d on detector %d is %d samples (should be %d)... skipping it.',evt,Detectors(i),length(DataEvt),ib));
            continue
        end
        
        DataEvt = DataEvt*SCALE; % now in Volts
        v = PassBand( DataEvt, ts, 50.e6, 100.e6 );
        vcor = spline(ts,v,tsover);  % Oversampling for better amplitude determination
        
        if ismember(nrun,runs2010) % Old soft
            fseek(ft,(evt-1)*ibufft*4,'bof');
            buff2=fread(ft,ibufft,'*uint32');
            trig = double(buff2(4));
        else
            trig = floor(ib/2);
        end
        
        % Time correction
        % Selection de la partie de la fonction d'onde autour du tps trig,
        % avec une largeur de 2*deltat.
        deltac= deltat*fover;
        trigc = trig*fover;
        v2 = double( vcor(max(1,trigc-deltac):min(trigc+deltac,length(vcor))) )';
        [ vpeakpeak, tmoy ] = FindVt( v2, trigc, deltac ); % calcul de vpeak et tmoy
        AmpMax(ind_tag(j),i) = vpeakpeak; % peak-peak amplitude in deltat window around trigger time
        TrigTime(ind_tag(j),i) = t - trig + tmoy/fover;   % Delay with first antenna in coinc
        
        %if DISPLAY
        if 0
            figure(1)
            plot(ts,v,'-b');
            hold on
            plot(tsover,vcor,'-r');
            disp(sprintf('Line %d - Coinc %d, Antenna %d: Amplitude = %3.2f V',ind_tag(j),CoincStruct.IdCoinc(ind_tag(j)),Detectors(i),vpeakpeak))
            pause
        end
        
    end
    disp(sprintf('Treatment completed for detector %d.',Detectors(i)));
end
for j=1:NbCoinc
    in = find(tag(j,:)==1);
    TrigTime(j,in) = TrigTime(j,in)-min(TrigTime(j,in));  %use 1st triggering detector as reference.
end

Struct.Coinc.Det.AmpMax = AmpMax;
Struct.Coinc.Det.TrigTime = TrigTime;

        
        
        
        
    
    
