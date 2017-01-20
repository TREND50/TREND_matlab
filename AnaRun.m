function AnaRun(nrun)
% Basic caras of runs
% OMH  02/06/11

SharedGlobals;
scrsz = get(0,'ScreenSize');

%% Load dst
dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
dst = load(dstname);

ncoincs = dst.Struct.Setup.TotalCoinc;
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
DetectorType = [dst.Struct.Setup.Det.isScint];
nDets = length(Detectors);
nScints = sum(DetectorType);
nAnts = nDets-nScints;
Events = [DetStruct.Evt];
CoincStruct = dst.Struct.Coinc;
tag = CoincStruct.Det.Tag;
L = CoincStruct.Mult;
id =  CoincStruct.Det.Id;
evt = CoincStruct.Det.Evt;
time = CoincStruct.Det.Time;
trig = CoincStruct.Det.TrigTime;
amp = CoincStruct.Det.AmpMax;
gain = CoincStruct.Det.Gain;
sat = [CoincStruct.Det.Sat];
stat = [CoincStruct.Det.Status];

InfoStruct = dst.Struct.Setup.InfosRun;
tstart = min(InfoStruct.TimeStart);
tstop = max(InfoStruct.TimeStop);
dur = (tstop-tstart)/60;
detTrigRate = [InfoStruct.TrigRate];
detCoincRateRaw = [InfoStruct.DetCoincRateRaw];
globalCoincRateRaw = [InfoStruct.GlobalCoincRateRaw];
detCoincRateCor = [InfoStruct.DetCoincRateQuickReject];
globalCoincRateCor = InfoStruct.GlobalCoincRateQuickReject;

disp(sprintf('Run %d - %3.1f mins - %d coincs',nrun,dur,ncoincs))
disp('***Analysis with triggers included in reconstructed coincidences only***')
for i=1:nDets
    if DetectorType(i)==0
        ampi = amp(:,i);
        gaini = gain(:,i);
        stdi = 1./gaini;
        sati = sat(:,i);
        satu = find(sati==1);
        tagi = tag(:,i);
        stati = stat(:,i);
        satandtrig = length(find(sati==1 & tagi==1));
        satandnotrig = length(find(sati==1 & tagi==0));
        in = find(tagi>0);
        t = time(:,i);
        t = (t-tstart)/FSAMPLING/60;  %time in minutes
        %
        % Amplitude
        if length(in)>0
            figure(1)
            set(1,'Name',sprintf('R%d Antenna%d - Amplitude',nrun,Detectors(i)),'NumberTitle', 'off')
            plot(t(in),ampi(in),'-')
            grid on
            hold on
            plot(t(satu),ampi(satu),'r.');
            xlim([0 dur]);
            xlabel('Time [min]', labelOpts{:})
            ylabel('Amplitude [LSB]', labelOpts{:})
            hold off
            %
            disp(sprintf('Antenna %d - %d triggers (%d saturated after treatment) - %d/%d in coincs - Av amplitude = %3.1f V - Av sigma = %3.1f V',Detectors(i),Events(i),sum(sati),sum(tagi),ncoincs,mean(ampi(in)),mean(stdi(in))))
            %disp(sprintf('%d saturated events in coincs',satandtrig))
            %
            %
            figure(2)
            set(2,'Name',sprintf('R%d Antenna%d - TriggerRate',nrun,Detectors(i)),'NumberTitle', 'off')
            subplot(2,1,1)
            plot(globalCoincRateRaw)
            hold on
            plot(globalCoincRateCor,'g')
            xlabel('Time [mn]')
            ylabel('Rate [Hz]')
            title('All coincidences')
            xlim([0 dur]);
            hold off
            subplot(2,1,2)
            plot(detTrigRate(:,i),'k')
            hold on
            plot(detCoincRateRaw(:,i))
            plot(detCoincRateCor(:,i),'g')
            title(sprintf('Coincidences with antenna %d',Detectors(i)))
            xlabel('Time [mn]')
            ylabel('Rate [Hz]')
            xlim([0 dur]);
            hold off
            pause
        end
        
        %
    end
end

