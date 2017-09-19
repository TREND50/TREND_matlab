function readNEvents(periodID)
% Simly checks nb of events per run in given period
% OMH 15/07/2017

SharedGlobals;

if periodID == 6
  runstart = 3562;
  runstop = 3563;
end

filename = [TEXT_PATH 'nevents.txt'];
fif = fopen( filename, 'a' );

for runid = runstart:runstop
  dstname = [DST_PATH sprintf(dst_filename,runid,1)];
  fid = fopen(dstname);
  if fid<0
      disp(sprintf('No dst %s',dstname))
      continue
  end
  fclose(fid);
  disp(sprintf('Loading dst %s...',dstname))
  dst = load(dstname);
    
    gcrr = dst.Struct.Setup.InfosRun.GlobalCoincRateRaw;
    rawcoincrate_mean = mean(gcrr(isfinite(gcrr)));  %[Hz]
    t = dst.Struct.Setup.InfosRun.TrigTime; %[s]
    if size(t,1) == 0
        disp 'No trigger rate measurment! Abort.'
        continue
    end
    tstart = t(1);
    tstop = t(end);
    durs = (tstop-tstart);  % [s]
    
    ncoincsraw = rawcoincrate_mean*durs;
    nevts = dst.Struct.Setup.TotalEvt;
    ncoincs = dst.Struct.Setup.TotalCoinc;
    
    fprintf( fif, '%6d ', runid );   % Event time after correction (in sample counts)
    fprintf( fif, '%8d ',  nevts );    % Antenna ID 
    fprintf( fif, '%3.1f ',  ncoincsraw );    % Event number
    fprintf( fif, '%8d ',  ncoincs );    % Event number % Wrong because only 1st DST is read!
	fprintf( fif, '\n' );
end

fclose(fif)