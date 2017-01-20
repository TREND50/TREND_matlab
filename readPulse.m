function [] = readPulse(nrun,antenna,evt)
% Simply reads and display signal of 
% event 'evt' on antenna 'antenna' in run 'nrun'.
% OMH 14/06/2011

SharedGlobals;
PATH = RAWDATA_PATH;
disp(sprintf('Looking for data in folder %s.',PATH))
filename = sprintf( 'R%05d_%d_data.bin' ,nrun, antenna );
fd = fopen( [PATH, filename ], 'r' );

%% Get data
if fd<0
    disp(sprintf('Could not find file %s.',filename))
else
    disp(sprintf('Opened file %s.',filename))
end
fseek( fd, ibuff*(evt-1)+1,'bof');
DataEvt = double( fread( fd, ibuff, 'uint8' ) );
t = [1:ibuff]/FSAMPLING;
tmu = t*1e6; %mus

%% Display
figure(1)
set(1,'Name',sprintf('Run %d Antenna %d Evt %d',nrun,antenna,evt),'NumberTitle','off')
plot(tmu,DataEvt,'k','LineWidth',2)
grid on
xlim([min(tmu) max(tmu)])
xlabel('Time [mus]', labelOpts{:})
ylabel('Amplitude [LSB]', labelOpts{:})
