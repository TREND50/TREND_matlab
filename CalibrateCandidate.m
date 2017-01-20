function [ calibSignals signals] = CalibrateCandidate( Struct )
% Perform calibration for shower candidates signals
% OMH 09/06/2011

RunSetup = Struct.Setup;
nrun = RunSetup.Run;
SharedGlobals;

Detectors=[RunSetup.Det.Name];
DetectorType=[RunSetup.Det.isScint];
CoincStruct = Struct.Coinc;
tag = CoincStruct.Det.Tag;
% Recons paras
IsShower = CoincStruct.IsShower;
CoincId = CoincStruct.IdCoinc;
ShowerId = CoincId(find(IsShower==1));
Lant = CoincStruct.MultAnt(find(IsShower==1));
L = CoincStruct.Mult(find(IsShower==1));
t = [ 1:ibuff ]/FSAMPLING; % seconds
tmu = t*1e6; % microsecs
fover = 10;
tover = [1 :1/fover: ibuff]/FSAMPLING;
tmuover = tover*1e6;
Calib = RunSetup.InfosRun.Calib;
freqPSDAll = Calib.F;

calibSignals = zeros(sum(L),3);
signals = zeros(sum(L),ibuff);

l=0;
for i=1:length(ShowerId)
    %disp(sprintf('Performing calibration for coinc %d',ShowerId(i)))
    ind_coinc = find(CoincId==ShowerId(i));
    ind_det = find(tag(ind_coinc,:)==1);
    ind_ant = find(tag(ind_coinc,:)==1 & DetectorType==0);
    Dets = CoincStruct.Det.Id(ind_coinc,ind_det);
    %Antennas = CoincStruct.Det.Id(ind_coinc,ind_ant);
    Events = CoincStruct.Det.Evt(ind_coinc,ind_det);
    GPSD = CoincStruct.Det.GainPSD{ind_coinc};
    a = 0;

    for j=1:length(Dets)
        l = l+1;                
        %disp(sprintf('Event %d Detector %d',Events(j),Dets(j)))
        % Get data
        fd = OpenFileData( nrun, Dets( j ) );
        if fd<0
            disp(sprintf('Could not find data for detector %d.',Dets(j)))
            fclose all
            continue
        end
        fseek( fd, ibuff*(Events(j)-1),'bof');
        DataEvt = double( fread( fd, ibuff, 'uint8' ) );
        DataEvt = DataEvt*SCALE; % Now in Volts
                
        
        %% Correct signal form  (method 2 using system response)
        if size(GPSD,1)==0
            calibSignals(l,1) = ShowerId(i);
            calibSignals(l,2) = Dets(j);
            signals(l,:) = DataEvt;
            continue % Skip calibration because no PSD is available
        end
        [freq,DataFFT,m] = FourierTrans(DataEvt,t,0);
        %rangeFFT = find(freq>=FREQMIN-5 & freq<=FREQMAX+5);
        %freq = freq(rangeFFT);
        %DataFFT = DataFFT(rangeFFT);
        a = a+1;
        GainP = GPSD(a,:);  % TO BE FIXED
        %rangePSD = find(freqPSD>=FREQMIN-5 & freqPSD<=FREQMAX+5);
        %freqPSD = freqPSD(rangePSD);
        %Gain = Gain(rangePSD);
        GainN = fliplr(GainP); % Symetrize
        Gain = [GainN GainP];
        GainLin = [10.^(Gain/20)]';
        if size(DataFFT,1)~= ibuff
            disp(sprintf('Error! Data size for event %d on detector %d is %d samples (should be %d)... skipping it.',Events(j),Dets(j),size(DataFFT,1),ibuff));
            fclose all
            continue
        end
        k1 = floor(ibuff/size(GainLin,1));
        GainLinInterp = interp(GainLin,k1);  % Interpolate
        out = find(GainLinInterp<100);
        DataFFTin = DataFFT./GainLinInterp; % complex FFT at input
        DataFFTin(out) = 1e-1000;
        modDataFFTin = abs(DataFFTin);        
        [aa,DataEvtin,m] = FourierTrans(DataFFTin,freq,1); % Re-invert axis. see FourierTrans
        %DataEvtin = ifft(fftshift(DataFFTin)); % Re-invert axis. see FourierTrans
        DataEvtin=PassBand(DataEvtin,t,FREQMIN,FREQMAX);
        DataEvtSpline = spline(t,DataEvtin,tover);  % Oversampling fo better amplitde determination
        
        % Compute max
        deltat = 10;
        deltac= deltat*fover;
        trig = floor(ibuff/2);
        trigc = trig*fover;
        DataZoom = double( DataEvtSpline(trigc-deltac:trigc+deltac) )';
        [ vpeakpeak, tmoy ] = FindVt( DataZoom, trigc, deltac ); % calcul de vpeak et tmoy
        calibSignals(l,1) = ShowerId(i);
        calibSignals(l,2) = Dets(j);
        calibSignals(l,3) = vpeakpeak; % peak-peak amplitude in deltat window around trigger time
        signals(l,:) = DataEvtin;
        rawAmp = max(DataEvt)-min(DataEvt);
        %disp(sprintf('Output amplitude = %3.1f V, Antenna amplitude = %3.1f muV, Gain = %3.1f dB',rawAmp,vpeakpeak*1e6,20*log10(rawAmp/vpeakpeak)))
        %k2 = size(DataEvt,1)/size(DataEvtin,1);
        %DataEvtinInterp = interp(DataEvtin,k2);
        %DataEvtinInterp=passband(DataEvtinInterp,t,50e6,100e6);
        
        % Display
        if DISPLAY
            freqP = freqPSDAll(j,:);
            freqN = -fliplr(freqP);
            freqPSD=[freqN freqP];
        
            figure(Dets(j))
            set(Dets(j),'Name',sprintf('Coinc %d - Event %d Antenna %d',ShowerId(i),Events(j),Dets(j)),'NumberTitle','off')
            subplot(3,2,1)
            plot(tmu,DataEvt)
            axis([0, max(tmu), 0, 2^NBITS*SCALE])
            grid on
            xlabel('Time [mus]')
            ylabel('Output signal [V]')
            subplot(3,2,2)
            semilogy(freq,abs(DataFFT))
            xlim([-FREQMAX FREQMAX])
            %xlim([0, FREQMAX])
            xlabel('Frequency [Hz]')
            ylabel('Output FFT [LSB]')
            subplot(3,2,3)
            
            plot(freqPSD,Gain)
            xlabel('Frequency [Hz]')
            ylabel('Gain [LSB/V]dB')
            xlim([-FREQMAX FREQMAX])
            %axis([0,FREQMAX, 0, 100])
            subplot(3,2,4)
            semilogy(freqPSD,GainLin)
            xlabel('Frequency [Hz]')
            ylabel('Gain [LSB/V]lin')
            %xlim([-FREQMAX FREQMAX])
            axis([0, FREQMAX, 100, 1e5])
            grid on
            subplot(3,2,5)
            semilogy(freq,modDataFFTin);
            xlabel('Frequency [Hz]')
            ylabel('Input FFT [V]')
            xlim([-FREQMAX FREQMAX])
            %xlim([0, FREQMAX])
            subplot(3,2,6)
            plot(tmu,DataEvtin,'LineWidth',1)
            hold on
            plot(tmuover,DataEvtSpline*1e6,'-m','LineWidth',1)
            grid on
            xlim([0, max(tmu)])
            xlabel('Time [mus]')
            ylabel('Input signal [muV]')
            pause
        end
        fclose(fd);
    end
end

