function [Struct]=ShowerLoader(Struct,calib)
% Load C++ shower reconstruction files and load them into dst
% 07/06/2011 OMH

SharedGlobals;


CoincId = Struct.Coinc.IdCoinc;

if ~exist('calib')
  calib = 0;
end

% Different case for subroutine used for normal delays or corrected delays
if ~isfield(Struct.Coinc,'DelayCorrRecons')    
    ref_struct=Struct.Coinc;
else
    ref_struct=Struct.Coinc.DelayCorrRecons;
end;

doShower = ref_struct.IsShower;

if calib==0
    ref_struct.ShowerRecons.RawRecons.XCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.RawRecons.YCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.RawRecons.ZCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.RawRecons.AxisAmp = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.RawRecons.Lambda = -1*ones(length(CoincId),1);
else
    ref_struct.ShowerRecons.CalRecons.XCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.CalRecons.YCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.CalRecons.ZCore = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.CalRecons.AxisAmp = -1*ones(length(CoincId),1);
    ref_struct.ShowerRecons.CalRecons.Lambda = -1*ones(length(CoincId),1);
end

%% Load recons files
nrun=Struct.Setup.Run;
filename = [TEXT_PATH sprintf( 'R%d_shower.txt', nrun)];
disp(sprintf('Opening file %s...',filename)); 
if fopen( filename )>0
    showerRes = load( filename );        
else
    disp(sprintf('File %s does not exist.',filename));
    if ~isfield(Struct.Coinc,'DelayCorrRecons')    
        Struct.Coinc=ref_struct;
    else
        Struct.Coinc.DelayCorrRecons=ref_struct;
    end;
    return
end  
disp 'Done.'
if size(showerRes,1)==0
    disp(sprintf('File %s empty.',filename));
    if ~isfield(Struct.Coinc,'DelayCorrRecons')    
        Struct.Coinc=ref_struct;
    else
        Struct.Coinc.DelayCorrRecons=ref_struct;
    end;
    return
end  
showerId = showerRes(:,1);
nCoincs = length(showerId);  % Number of reconstructed showers
xc = showerRes(:,3);
yc = showerRes(:,4);
zc = showerRes(:,5);
amp0 = showerRes(:,6);
lambda0 = showerRes(:,7);


%% Write to structure
disp('Now writting shower reconstruction results to DST...')

for i = 1:nCoincs
    ind = find(CoincId==showerId(i));  % Index of shower in reconstruction structure
    if size(ind,2)==0
        disp(sprintf('Error in shower reconstruction for coinc %d',showerId(i)))
        continue
    end
    if doShower(ind) == 1
        if calib == 0
            ref_struct.ShowerRecons.RawRecons.XCore(ind) = xc(i);
            ref_struct.ShowerRecons.RawRecons.YCore(ind) = yc(i);
            ref_struct.ShowerRecons.RawRecons.ZCore(ind) = zc(i);
            ref_struct.ShowerRecons.RawRecons.AxisAmp(ind) = amp0(i);
            ref_struct.ShowerRecons.RawRecons.Lambda(ind) = lambda0(i);
        else
            ref_struct.ShowerRecons.CalRecons.XCore(ind) = xc(i);
            ref_struct.ShowerRecons.CalRecons.YCore(ind) = yc(i);
            ref_struct.ShowerRecons.CalRecons.ZCore(ind) = zc(i);
            ref_struct.ShowerRecons.CalRecons.AxisAmp(ind) = amp0(i);
            ref_struct.ShowerRecons.CalRecons.Lambda(ind) = lambda0(i);
        end
    else
        disp(sprintf('Error in shower reconstruction for coinc %d',showerId(i)))
    end
end

if ~isfield(Struct.Coinc,'DelayCorrRecons')    
    Struct.Coinc=ref_struct;
else
    Struct.Coinc.DelayCorrRecons=ref_struct;
end;