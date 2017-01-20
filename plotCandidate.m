function [  ] = plotCandidate( nrun,id, ddst )
% Displays signals for a shower candidate
% OMH 21/06/2011

SharedGlobals;
scrsz = get( 0, 'ScreenSize' );

%% Plot options
%%===
labelOpts = { 'FontSize', 14 };
axisOpts  = { 'fontSize', 14, 'YScale', 'linear' };

%% Load dst
if ~exist('ddst')
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    ddst = load(dstname);
end
Struct = ddst.Struct;
CoincStruct = Struct.Coinc;
DetStruct = Struct.Setup.Det;
Detectors = [DetStruct.Name];

% ShowerPars
CoincId = CoincStruct.IdCoinc;
tag = CoincStruct.Det.Tag;
Evt = CoincStruct.Det.Evt;
IsShower = CoincStruct.IsShower;
signals = CoincStruct.ShowerRecons.ShowerSignals;
L = CoincStruct.Mult;
ind = find(CoincId==id);
if size(ind,1)==0
    disp(sprintf('Coinc %d not reconstructed. Aborting.',id))
    return; 
elseif IsShower(ind)==0
    disp(sprintf('Coinc %d is not a shower candidate. Aborting.',id))
    return;
end
indant = find(tag(ind,:)==1);

t = [1:ibuff]/FSAMPLING;
tmu = t*1e6;
for i = 1:L(ind)
    rawEvt = signals{ind,indant(i)}*1e6;  % Amplitude in muV
    if length(rawEvt)~=ibuff
        disp(sprintf('Error! Data size in coinc %d for detector %d is %d samples (should be %d)... skipping it.',id,Detectors(indant(i)),length(rawEvt),ibuff));
        continue
    end
    figure(1)
    set(1,'Name',sprintf('Run %d Coinc %d ',nrun,id),'NumberTitle','off')
    if L(ind)==4
      nlin = 2;
      ncol = 2;
    elseif L(ind)<=6
        nlin = 3;
        ncol = 2;
    elseif L(ind)<=9
        nlin=3;
        ncol = 3;
    elseif L(ind)<=16
        nlin=4;
        ncol = 4;
    else
        nlin = 5;
        ncol = 5;
    end
    subplot(nlin,ncol,i)
    plot(tmu,rawEvt,'k')
    grid on
    xlim([0 max(tmu)])
    xlabel('Time [mus]')
    ylabel('Amp [muV]')
    title(sprintf('Antenna %d - Event %d',Detectors(indant(i)),Evt(ind,indant(i))),labelOpts{:})
end

