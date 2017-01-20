function AnaBckg( periodId, nrun )
% Determine rate of events triggering antennas in [antennas] vector.
% To be used for GRAND bckgrd rate estimation.
% OM 17/04/2016

if periodId>0
    switch (periodId)
      case 1
        nruns = 2538:2585;
      case 2
        nruns = 2685:2890;
      case 3
        nruns = 3000:3086;
      case 4
        nruns = 3157:3256;        
      case 5
        nruns = 3336:3371;   
      case 6
        nruns = 3562:3733;
      case 7
        nruns = 3835:3999;
      case 8
        nruns = 4183:4389;
      case 9
        nruns = 4444:5014;
      case 10
%        nruns = 5070:5913;
        nruns = 5071:5913;
    end
else
    nruns =nrun;  % Provide run list by hand
end
logname = sprintf('coinctable_GRAND_%d.txt',periodId);
durtot = 0;
dur = 0;
tlive = 0;
tlivetot = 0;
ncoinctot = 0; 
for i=1:length(nruns)
    irun = nruns(i)
    if irun<3561
        disp('Error: log files format before R3561 does not allow to compute live time.')
        %return
    end

    if ~exist('dis')
        dis = 1;
    end
    SharedGlobals;
    if periodId == 9
	  antennas = [101 111 119 128 150 158]  
    else
	  antennas = [101 111 119 128 149 158]  
    end
    
    
    if irun>=3561
        %% Compute live time
        [tnotrig tlate tsat tdead tmin tmax] = getLogData(irun,antennas);
        if size(tlate,2)<length(antennas)
            disp 'At least one antenna log file is missing. Discard this run.'
            continue
        end
        %
        trun = round(tmin):round(tmax);
        dur = length(trun)/3600;
        for j = 1:length(antennas)
           trun = setdiff(trun,tlate{j});
           trun = setdiff(trun,tsat{j});
           trun = setdiff(trun,tnotrig{j});
           trun = setdiff(trun,tdead{j});
           %length(trun)
        end
        tlive = length(trun)/3600;
        disp(sprintf('R%d: duration = %3.1f hours, live time with all antennas = %3.1f hours.',irun,dur,tlive))
        if dur>100 %more than 100hours
           disp 'Error with log file time'
           break
        end
    end
    
    %% Get nb of coincs
    ncoincs = 0;
    for k = 1:10
        dstname = [DST_PATH sprintf('dst%d_%d.mat',irun,k)];
        disp(sprintf('Loading dst %s...',dstname))
        if fopen(dstname)<0
            disp 'No dst available. Abort.'
            break
        else
            dst = load(dstname);
            disp 'Done.'
        end
        dets = [dst.Struct.Setup.Det.Name];
        [c iants] = intersect(dets,antennas);
        tag = dst.Struct.Coinc.Det.Tag(:,iants);  % 
        %tag
        %sum(tag,2)
        allin = find(sum(tag,2)==length(antennas)); % All target antennas are in this coinc
        %
        ut = dst.Struct.Coinc.Det.UnixTime(allin,iants);  % Get Unix Time for these coincs
        ut = round(mean(ut,2));
        if periodId>=6
            allgood = setdiff(ut,trun);  % Check if in list of 'good' seconds
            if length(allgood)>0
                disp 'Error!'
                allgood
            end
        else % Older periods: compute duration from DST time info.
          tlive = 0;  
          tstart = median(dst.Struct.Setup.InfosRun.TimeStart);
          tstop = median(dst.Struct.Setup.InfosRun.TimeStop);
          dur = (tstop-tstart)/3600;  % [h]
        end
        %
        id = dst.Struct.Coinc.IdCoinc(allin);
        l = dst.Struct.Coinc.MultAnt(allin);
        thp = dst.Struct.Coinc.PlanRecons.Radio.Theta(allin);
        phip = dst.Struct.Coinc.PlanRecons.Radio.Phi(allin);
        chi2p = dst.Struct.Coinc.PlanRecons.Radio.Chi2Delay(allin);
        slopep = dst.Struct.Coinc.PlanRecons.Radio.SlopeDelay(allin);
        xs = dst.Struct.Coinc.SphRecons.X0(allin);
        ys = dst.Struct.Coinc.SphRecons.Y0(allin);
        zs = dst.Struct.Coinc.SphRecons.Z0(allin);
        chi2s = dst.Struct.Coinc.SphRecons.Chi2Delay(allin);
        slopes = dst.Struct.Coinc.SphRecons.SlopeDelay(allin);
        
        evtname = sprintf('coinctable_GRAND_%d_details.txt',periodId);
        fid = fopen( evtname, 'a+' );   
        for ii = 1:length(allin)
          fprintf( fid, '%6d %6d %12.3f %6d %3.1f %3.1f %5.1e %3.2f %3.5e %3.5e %3.5e %5.1e %3.2f \n',  irun, id(ii), ut(ii), l(ii), thp(ii), phip(ii), chi2p(ii), slopep(ii), xs(ii), ys(ii), zs(ii), chi2s(ii), slopes(ii));     
        end
        fclose(fid);
        ncoincs = ncoincs+length(ut);
    end
    disp(sprintf('%d coincs with these antennas in this run.',ncoincs))
    %
    durtot = durtot + dur;
    tlivetot = tlivetot+tlive;
    ncoinctot = ncoinctot+ncoincs;
    
    fid = fopen( logname, 'a+' );        
    fprintf( fid, '%6d %3.1f %3.1f %6d\n',  irun, dur, tlive, ncoincs);     
    fclose(fid);

  end
  durtot
  tlivetot
  ncoinctot
  rate = ncoinctot/(tlivetot*3600)
end

function [tnotrig tlate tsat tdead tmin tmax] = getLogData(nrun,antennas)
    SharedGlobals;
    tdead = {};
    tlate = {};
    tnotrig = {};
    tsat = {};
    tmin = 1e20;
    tmax = 0;
    for i = 1:length(antennas)
        antenna = antennas(i);
        filename = sprintf( 'R%06d_A%04d_log.txt', nrun, antenna ); 
        filename = [LOG_PATH, filename];
        if fopen(filename)<0
            disp(sprintf('Could not find file %s',filename))
            continue
        end
        logfile=load(filename);
        tloop = logfile(:,1); % Time at begining of reading loop
        irqs =  logfile(:,7);
        irqe =  logfile(:,8);
        nspike = logfile(:,9);
        ntrig = logfile(:,10);
        lsb = logfile(:,11);
        %
        tmin = min(tmin,min(tloop));
        tmax = max(tmax,max(tloop));
        %
        check = irqe-irqs;
        ilate = find(check>=1);
        tlate{i} = round(tloop(ilate));
        %
        isat = find(nspike == 256);
        tsat{i} = round(tloop(isat));
        %
        idead = find(lsb<1);
        tdead{i} = round(tloop(idead));
        %
        win = 180;  % Window to check that coincs exist [s]
        inotrig = [];
        for j = 1:win:length(tloop)-win;
           if (sum(ntrig(j:j+win-1)))==0 
               inotrig(end+1:end+win) = j:j+win-1;
           end
        end
        tnotrig{i} = round(tloop(inotrig));
        %
%         length(inotrig)
%         length(idead)
%         length(ilate)
        fclose all;
    end

end
