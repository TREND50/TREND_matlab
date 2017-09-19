function [] = computeTrigRate(runid)
% Loads RUNID dst and reads trigger rate
% Inserts in in total_trigrate and writes to trigrate.mat
% OMH 27/05/2013

SharedGlobals;
DISPLAY = 1;

%% Load structure
filename = ([MONITOR_PATH 'trigrate.mat']);
fid = fopen(filename);
if fid>0
    disp(sprintf('Loading %s',filename))
    t = load(filename);
    unixs = t.unixs;
    rate = t.coincrate;
    trate = t.t0rate;
    allruns = t.allruns;
    disp 'Done.'
    fclose(fid);
else
    unixs = [];
    rate = [];
    trate = zeros(50,1);
    allruns = {};
end

runs = [];
nruns = size(allruns,2);
for i =1:nruns
  runs(end+1) = allruns{i}.id; 
end    

%% Deal with run runid
if sum(ismember(runs,runid))>0
    disp(sprintf('Run %d already in %s. Nothing to do. ',runid,filename))
else
    disp(sprintf('Now processing run %d... ',runid))
    %% Read dst
    dstname = [DST_PATH sprintf(dst_filename,runid,1)];  % Only needed to read 1st dst, as we load Setup substruct only
    fid2 = fopen(dstname);
    if fid2<0
        disp(sprintf('dst %s not found.',dstname))
    else
        fclose(fid2);
        disp(sprintf('Loading dst %s...',dstname))
        dst = load(dstname);
        disp 'Done.'
        t = dst.Struct.Setup.InfosRun.TrigTime; %[s]
        if size(t,1) == 0
            disp 'No trigger rate measurment! Abort.'
            return
        end
        tstart = t(1);
        tstop = t(end);
        durs = (tstop-tstart);  % [s]
        dur = durs/60; %[mn]
        durh = dur/60; %[h]        
        info = dst.Struct.Setup.InfosRun;
        if ~isfield(info,'MetaTotalCoinc')
            disp 'No field MetaTotalCoinc' 
            nt1 = 0;
            ntf = 0;
        else
            nt1 = sum(info.MetaTotalCoinc);
            ntf = sum(info.TotalCoincFiltered2);
        end
        
        coincrate = info.GlobalCoincRateRaw;
        nt = size(coincrate,1);
        if nt>10000
            disp(sprintf('Run %d: %3.1f h and %d points! Problem here. Abort.',runid,durh,nt))
            return
        end
        step = round(dur/nt); % [mn]
        if step<1
            disp(sprintf('Error in computing trigger rate step: step = round(%3.1f/%d)=%d. Abort.',dur,nt,step))
            return
        end
        steph = step*60;
        disp(sprintf('Run %d: %3.1f h, %d points (%3.2f / mn)',runid,durh,nt,step))
        detrate = info.DetCoincRateRaw;
        det0rate = info.TrigRate;
        nevts = [dst.Struct.Setup.Det.Evt];
        dets = [dst.Struct.Setup.Det.Name];    
        ndets = length(dets);

        %% Now fill in vectors
        i = 1;
        % Event number
        % Trig rate
        th = zeros(1,floor(nt/steph));
        detratec = zeros(1,length(dets));
        coincrateh = zeros(1,floor(nt/steph));     
        t0rateh = zeros(length(ALLDETS),floor(nt/steph));     
        
        while i<=max(1,nt/steph)
            indbeg = (i-1)*steph+1;
            indend = min(i*steph,nt);
            th(i) = tstart+(i-1)*3600; % [s]
            stepint = 60/step; %time interval between 2 rate measurment [s]
            coincratei = coincrate(indbeg:indend);
            coincrateh(i) = sum(coincratei(isfinite(coincratei))*stepint)/3600; 
            for j = 1:ndets
              indd = find(ALLDETS == dets(j));  
              t0ratei = det0rate(indbeg:indend,j);
              t0rateh(indd,i) = sum(t0ratei(isfinite(t0ratei))*stepint)/3600;
            end
            % AvTrigRate over time slot [Hz] = Sum of events in this time slot / time slot duration [s]
            disp(sprintf('Slice %d: mean coinc rate = %3.1f Hz',i,coincrateh(i)))
            % Detector data
            detnbt1 = nevts; % Nb of T1 events for each antenna
            % Average trigger rate for events in coincs each antenna on this run
            detratec = sum(detrate(:,:)*stepint,1)/durs;  
            i = i+1;
        end

        %% Now complete structure
        run.id = runid;
        run.dets = dets;
        run.dur = durh;
        run.nt1 = nt1;
        run.ntf = ntf;
        run.detcrate = detratec;
        run.detnbt1 = detnbt1;      
        allruns{end+1} = run;
        %
        if fid>0   %trigrate.mat exists
          [a ind] = min(abs(unixs-th(1)));
          unixs = [unixs(1:ind-1) th unixs(ind:end)];
          rate = [rate(1:ind-1) coincrateh rate(ind:end)];
          for j = 1:length(ALLDETS)
              newcol = [trate(j,1:ind-1) t0rateh(j,:) trate(j,ind:end)];
              newtrate(j,1:length(newcol)) = newcol;            
          end
          trate = newtrate;
        else   %trigrate.mat does not exist
          unixs = th;
          rate = coincrateh;
          for j = 1:ndets
            indd = find(ALLDETS==dets(j));  
            if size(indd,2)>0
                tamere=t0rateh(j,:);
                trate(indd,1:length(tamere)) = tamere;
            end
          end
        end
        coincrate = rate;
        t0rate = trate;
        save(filename,'allruns','unixs','coincrate','t0rate')
    end
end

