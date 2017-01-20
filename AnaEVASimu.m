function [trig ispres theta phi pol] = AnaEVASimu(jobid,E,xCore,yCore)
% Plot results of EVA & NEC2 simulations
% Add noise & compute trigger
% Data saved to *.bin files
% Analysis results saved to mat files.

% OMH 24/11/2012 EVA
% + 22/02/2013 NEC2
% + 16/12/2013 new EVA + adapted to TREND-50
% + 24/01/2014 adapted to new data structure for simulation
% + 24/03/2014 adapted to MassProd data: now no more loop on CorePos.

  SharedGlobals;
  DISPLAY = 0;
  N = 8;  % Threshold
  
  col = {'k','b','r','k','b','r'};
  tit = {'E_{EW}','E_{NS}', 'E_{vert}','EW antenna voltage','NS antenna voltage','Vert antenna voltage'};
  f = {'efield.txt','ew.txt','ns.txt','vert.txt'};
  
  %% Load coordinates
  coord = load('coord_antennas_all.txt');
  idant = coord(:,1);
  x = coord(:,2);
  y = coord(:,3);
  z = coord(:,4);
  id = 101:158;  % TREND-50 antennas
  
  amp = zeros(1,6);
  pol = zeros(1,7);
  trig = zeros(1,3);
  ispres = zeros(1,3);
  theta = 0;
  phi = 0;

  %% Load DST
  if E == '1e17'
    nrun = 1017;
  elseif E == '5e17'
    nrun = 5017;
  elseif E == '8e16'
    nrun = 8016;
  elseif E == '5e16'
    nrun = 5016;
  elseif E == '2e17'
    nrun = 2017;
  end 
  filename =sprintf(dst_filename,nrun,1);
  if fopen(filename)>0
      load(filename);
      dstDet = Struct.Setup.Det;
      dstCoinc = Struct.Coinc;
      dstSimu = Struct.Simu;       
      dstNbOfCoinc = Struct.Setup.TotalCoinc;
  else
      dstDet = {};
      dstDet.Name = [];
      dstNbOfCoinc = 0;
      dstCoinc = {};
      dstCoinc.Det = {};
      dstCoinc.Det.Id = [];
      dstCoinc.Det.Tag = [];
      dstCoinc.Det.Evt = [];
      dstCoinc.Det.Cal = [];
      dstCoinc.Det.TrigTime = [];
      dstCoinc.Det.TrigCor = [];
      dstCoinc.Det.AmpMax = [];
      dstCoinc.Det.Sigma = [];

      dstCoinc.Mult = [];
      dstCoinc.MultAnt = [];
      dstCoinc.IdCoinc = [];
      dstCoinc.IsShower = [];
      dstSimu = {};
  end
  
  necfolder = sprintf('/%d*%d_%d*',jobid,xCore,yCore);
  necfolder = [SIMU_PATH E '/voltages/' necfolder];  
  %fl = ls(necfolder);
  %clear fl;
  dd = dir(necfolder);
  if size(dd,1)==0 
     disp(sprintf('No folder %s. Abort.',necfolder))
      %disp(sprintf('Existing folders with shower id = %d:',jobid))
      %ls(sprintf('%s%d*',SIMU_PATH,jobid))
     fclose all;
     clear necfolder;
     clear dd;
     return
  end
  fl = dd(1).name;

  eta = zeros(1,length(id));
  beta = zeros(1,length(id));
  eta_v = zeros(1,length(id));
  beta_v = zeros(1,length(id));
  
  %
  truename= strread(fl(:),'%s','delimiter',' ');
  clear fl;
  clear dd;
  fclose all;
  folder = sprintf('%s/%s/voltages/%s/',SIMU_PATH,E,truename{1});
  [ast theta phi_eva xCore zst] = strread(folder,'%s%d%d%d%s','delimiter','_');
  phi = mod(phi_eva-90,360);
  zCore = REFALT;

  tic
  disp(sprintf('Zen = %d deg, Az = %d deg, core position = (%dm, %dm)',theta,phi,xCore,yCore))
  maxi = zeros(length(id),7);
  txtTable = zeros(1,12);
  %ok = 0;

  %% Loop on antennas  
  for i = 1:length(id)
      ind = find(idant==id(i));
      %disp(sprintf('Antenna %d: (%3.0f,%3.0f,%3.0f)',idant(ind),x(ind),y(ind),z(ind)))

      %% Loop on channels
      for j=1:3
        % Loop on fiel=ds, with:
        % 1: Efield
        % 2: EW
        % 3: NS
        % 4: Vert

%             if ~(j==1 || j== channel+1)  % Consider only one channel when studying trigger
%                 continue
%             end

        %% Load files
        filename = sprintf('a%d_%s',id(i),f{j});
        filename = [folder filename];
        if fopen(filename)<0
          %disp(sprintf('No file %s. Skip.',filename))
          continue
        end
        % Check file format
        if j>1
           [s,m] = unix(sprintf('wc -l %s',filename));
           nl = str2num(m(1));
           if nl ~= 5
             disp(sprintf('%s: corrupted text file! Skip.',filename))
             continue
           end
        end

        m = load(filename);  
        if size(m,1) == 0
          %disp(sprintf('File %s is empty. Skip.',filename))
          continue
        end
        %           
        if j == 1  % Efield
          t = m(:,1)*1e-6; % s
          t_off = t(1)*FSAMPLING; %sample units
          xe = m(:,2); % muV/m
          ye = m(:,3); % muV/m
          ze = m(:,4); % muV/m
          %
          mx = abs(max(xe)-min(xe));
          my = abs(max(ye)-min(ye));
          mz = abs(max(ze)-min(ze));   
          %ok = ok+1;
          maxi(end+1,1) = id(i);
          maxi(end,2) = mx;
          maxi(end,3) = my;
          maxi(end,4) = mz;
          
          [m indm] = max(abs(xe));
          pol(1) = pol(1)+sign(xe(indm));
          [m indm] = max(abs(ye));
          pol(2) = pol(2)+sign(ye(indm));
          [m indm] = max(abs(ze));
          pol(3) = pol(3)+sign(ze(indm));
          
        else % Voltages
          ispres(j-1) = ispres(j-1)+1;    % data is present
          
          m = m';         
          if size(m,2)~=5 | size(m,1)<100
            disp 'Corrupted text file! Skip.'
            continue
          end
          t = m(:,1); % s
          vs = m(:,2); % muV
          v_r = m(:,3); % muV
          v_th = m(:,4); % muV
          v_ph = m(:,5); % muV
          v(j-1,:) = vs;
      %    if size(vs,1)<100
            
          maxv(j-1) = max(v(j-1,1:end-200))-min(v(j-1,1:end-200));
          [m it ] = max(abs(v(j-1,1:end-200)));
          tmax0 = t(it)*FSAMPLING;  % Samples
          maxi(end,j+3) = maxv(j-1);
          
          % Polarity
          [m indm] = max(abs(vs));
          pol(j+2) = pol(j+2)+sign(vs(indm));    % data is present

          % Shape signal
          offset = 500;% Offset by 500 points
          vs = [zeros(1,offset) vs(1:end-offset)']';  
          [vf tf] = PassBand(vs,t,FREQMIN,FREQMAX);  % Filter signal in system bandpass (50-100MHz)
          tf = tf*1e9;  %ns
          [td vd] = computeDigitizedSignal(tf,vf);  % Digitize signal at sampling frequency

          % Add Noise
          % Get experimental model
          % Load calib DST & noise
          c = load([CAL_PATH 'calib3135_1.mat']);
          ant_cal = c.Struct.Setup.Antennas;
          goodants = [108,111,113,123,126,149];  % Antennas with reliable results in 3129-3130 calib runs + noise level OK
          rind = ceil(rand(1)*length(goodants)); 
          acal = goodants(rind);
          %acal = 126;
          iind = find(ant_cal == acal);
          g = c.Struct.Calib.GainLin(iind)/SCALE; % V@antenna -> LSB@DAQ
          % Load backrground run (noise)
          fd = OpenFileBackground(3135,acal);
          data=fread(fd);
          fclose(fd);
          size_event = length(vd);
          start = round((size(data,1)-size_event)*rand(1,1));
          noise = data(start+1:start+size_event);
          th = N*std(noise);
          vdaq = vd*g*1e-6;
          vfin = vdaq+noise;

          % Saturation ADC
          vfin(vfin>255) = 255;
          vfin(vfin<0) = 0;   
          
          %% Trigger
          [mf it] = max(abs(vfin-mean(vfin)));
          if mf>th
              %disp(sprintf('Channel %d of antenna %d triggered!',j-1,idant(ind)))
              trig(j-1) = trig(j-1)+1;
              % Polarity
              vpol = vfin-mean(vfin);
              [m indm] = max(abs(vpol));
              pol(j+4) = pol(j+4)+sign(vpol(indm));    % data is present
              
              % Write to file (for reconstruction)
              if j==2  %only for EW
                  tmax = tf(it)/5; % Units of samples
                  if sum(ismember(txtTable(:,2),idant(ind)))==0
                        txtTable(end+1,:) = [0 idant(ind) 1 1 t_off+tmax t_off+tmax 0 0 0 0 0 0];
                  end
             
                  % Write to dstDet
                  indant = find([dstDet.Name]==idant(ind));
                  if size(indant,2) == 0  % Detector not yet in DST
                        if size(dstDet(1).Name,1) == 0  % Empty DST
                            dstDet(end).Name = idant(ind);
                        else
                            dstDet(end+1).Name = idant(ind);
                        end
                        dstDet(end).Evt = 1;
                        dstDet(end).X = x(ind);
                        dstDet(end).Y = y(ind);
                        dstDet(end).Z = z(ind);
                        dstDet(end).PodNb = idant(ind);                    
                        dstDet(end).PodX=0;
                        dstDet(end).PodY=0;
                        dstDet(end).PodZ=0;
                        dstDet(end).Cable=0;
                        dstDet(end).Delay=0;
                        dstDet(end).DelayCorr=0;
                        dstDet(end).isScint=0;     
                        dstDet(end).isLog=0;
                        %
                        if trig(1) == 1  % First triggered antenna for this coinc, add one line to Coinc matrices
                            dstCoinc.Det.Tag(end+1,end+1) = 1;
                            dstCoinc.Det.Evt(end+1,end+1) = 1;
                            dstCoinc.Det.Cal(end+1,end+1) = acal;                        
                            dstCoinc.Det.Id(end+1,end+1) = idant(ind);                            
                            dstCoinc.Det.TrigTime(end+1,end+1) = t_off+tmax;
                            dstCoinc.Det.TrigCor(end+1,end+1) = t_off+tmax;       
                            dstCoinc.Det.AmpMax(end+1,end+1) = mf;      
                            dstCoinc.Det.Sigma(end+1,end+1) = std(noise);
                        else
                            dstCoinc.Det.Tag(end,end+1) = 1;
                            dstCoinc.Det.Evt(end,end+1) = 1;    
                            dstCoinc.Det.Cal(end,end+1) = acal;                        
                            dstCoinc.Det.Id(end,end+1) = idant(ind);                            
                            dstCoinc.Det.TrigTime(end,end+1) = t_off+tmax;
                            dstCoinc.Det.TrigCor(end,end+1) = t_off+tmax;     
                            dstCoinc.Det.AmpMax(end,end+1) = mf;                             
                            dstCoinc.Det.Sigma(end,end+1) = std(noise);
                        end
                  else
                        dstDet(indant).Evt = dstDet(indant).Evt+1;
                        if trig(1) == 1  % First triggered antenna for this coinc, add one line to Coinc matrices
                            dstCoinc.Det.Tag(end+1,indant) = 1;
                            dstCoinc.Det.Evt(end+1,indant) = dstDet(indant).Evt;
                            dstCoinc.Det.Cal(end+1,indant) = acal;                        
                            dstCoinc.Det.Id(end+1,indant) = idant(ind);                            
                            dstCoinc.Det.TrigTime(end+1,indant) = t_off+tmax;
                            dstCoinc.Det.TrigCor(end+1,indant) = t_off+tmax;                            
                            dstCoinc.Det.AmpMax(end+1,indant) = mf;   
                            dstCoinc.Det.Sigma(end+1,indant) = std(noise);
                        else
                            dstCoinc.Det.Tag(end,indant) = 1;
                            dstCoinc.Det.Evt(end,indant) = dstDet(indant).Evt;
                            dstCoinc.Det.Cal(end,indant) = acal;
                            dstCoinc.Det.Id(end,indant) = idant(ind);                            
                            dstCoinc.Det.TrigTime(end,indant) = t_off+tmax;
                            dstCoinc.Det.TrigCor(end,indant) = t_off+tmax; 
                            dstCoinc.Det.AmpMax(end,indant) = mf;
                            dstCoinc.Det.Sigma(end,indant) = std(noise);
                        end
                  end
                  names = [dstDet.Name];
                  [a iord] = sort(names);
                  dstDet = dstDet(iord);    
                  dstCoinc.Det.Tag = dstCoinc.Det.Tag(:,iord);
                  dstCoinc.Det.Evt = dstCoinc.Det.Evt(:,iord);
                  dstCoinc.Det.Cal = dstCoinc.Det.Cal(:,iord);
                  dstCoinc.Det.Id = dstCoinc.Det.Id(:,iord);
                  dstCoinc.Det.TrigTime = dstCoinc.Det.TrigTime(:,iord);
                  dstCoinc.Det.TrigCor = dstCoinc.Det.TrigCor(:,iord);
                  dstCoinc.Det.AmpMax = dstCoinc.Det.AmpMax(:,iord);
                  dstCoinc.Det.Sigma = dstCoinc.Det.Sigma(:,iord);
              end
              
%             % Write data to file
%            filename = sprintf( 'S%s_A%04d_data.bin', fiid, idant(ind) ); 
%            fdw = fopen( [ folder, filename ], 'w+' );
%            if length(vfin)<ibuff
%              start2 = round((size(data,1)-ibuff)*rand(1,1));
%              vw = [data(start2:start2+ibuff-length(vfin))' vfin'];
%            end
%            fwrite(fdw,vw,'*uint8');
          
          else
              %disp(sprintf('Channel %d of antenna %d did not trigger.',j-1,idant(ind)))
          end
        end
        t = t-t(1);
        tns = t*1e9;  %ns

        %% Plot
        if DISPLAY == 2 
              if j == 1
                  figure(i*100+1)
%                  plot(tns,xe,'k','LineWidth',2)
                  hold on
                  grid on
%                   plot(tns,ye,'b','LineWidth',2)
%                   plot(tns,ze,'r','LineWidth',2)
                  plot(tns,xe+ye+ze,'k','LineWidth',2)
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Efield [muV/m]', labelOpts{:})
                  title(sprintf('Antenna %d - Efield',id(i)),labelOpts{:})
                  hold off
              else
                  figure(i*100+j)
                  %plot(tns,v_r,'k--','LineWidth',2)
                  hold on
                  grid on
                  %plot(tns,v_th,'b--','LineWidth',2)
                  %plot(tns,v_ph,'r--','LineWidth',2)
                  %plot(tns,v(j-1,:),'g--','LineWidth',2)
                  %plot(tns,v_r+v_th+v_ph,'y','LineWidth',2) %Check
                  %
                  plot(tns,vf,'k','LineWidth',1)
                  plot(td,vd,'mo-','LineWidth',1,'MarkerFaceColor','w')
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Voltage [muV]', labelOpts{:})
                  title(sprintf('Antenna %d - %s',id(i),tit{j+2}),labelOpts{:})
                  hold off
                  %
                  figure(i*1000+j)
                  subplot(2,1,2)
                  plot(td,vfin,'k','LineWidth',2)
                  grid on
                  mn = mean(noise);
                  th_up = line([min(td) max(td)],[mn+th mn+th]);               
                  th_down = line([min(td) max(td)],[mn-th mn-th]);
                  set(th_up,'LineStyle','--','Color','k')
                  set(th_down,'LineStyle','--','Color','k')        
                  ylim([min(vfin)*0.9 max(vfin)*1.1])
                  ylim([0 255])
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('DAQ signal+noise [LSB]', labelOpts{:})
                  title(sprintf('Antenna %d - %s',id(i),tit{j+2}),labelOpts{:})
                  subplot(2,1,1)
                  plot(td,noise,'r','LineWidth',2)
                  hold on
                  plot(td,vdaq+mean(noise),'b--','LineWidth',1)
                  ylim([min(vfin)*0.9 max(vfin)*1.1])
                  ylim([0 255])
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('DAQ signal [LSB]', labelOpts{:})
                  grid on
              end
              if j==4
                figure(i*100+1)
                pause
                close all  
              end
        end

        %% ETA & Beta computation
        if j==1
            plane=sqrt(xe.*xe+ye.*ye);
            mplane = abs(max(plane)-min(plane));
            eta(i) = atand(my/mx);
            beta(i) = atand(mplane/mz);
        elseif j==4
            plane=sqrt(v(1,:).*v(1,:)+v(2,:).*v(2,:));
            mplane = abs(max(plane)-min(plane));
            eta_v(i) = atand(maxv(2)/maxv(1));
            beta_v(i) = atand(mplane/maxv(3));
        end  

      end   % loop on channels         
      fclose all;
      clear v;
      clear v_r;
      clear v_th;
      clear v_ph;

  end  % loop on antennas

  if DISPLAY == 2
      pause
  else
    close all
  end

  trig
  pol

  %% Build output elements if trigger
  % Find triggering channels
  if trig(1)>=5  %EW channel triggered
        % Coinctable
        txtTable(1,:) = [];
        txtTable = sortrows(txtTable,5);
        txtTable(:,5) = txtTable(:,5)-txtTable(1,5);
        txtTable(:,6) = txtTable(:,6)-txtTable(1,6);
        %
        %
        % Build coinctable file
        if E=='1e17'
          filename = [TEXT_PATH sprintf( 'R1017_coinctable.txt')];
        elseif E=='5e17'
          filename = [TEXT_PATH sprintf( 'R5017_coinctable.txt')]; 
        elseif E=='8e16'
          filename = [TEXT_PATH sprintf( 'R8016_coinctable.txt')];
        elseif E=='5e16'
          filename = [TEXT_PATH sprintf( 'R5016_coinctable.txt')];
        elseif E=='2e17'
          filename = [TEXT_PATH sprintf( 'R2017_coinctable.txt')];
        end

        fid = fopen( filename, 'a+' );
        c = load(filename);
        if size(c,1)>0
          aid = c(:,2);
          cid = c(end,4);            
        else
          aid = 0;
          cid = 0;
        end
        for l = 1:size( txtTable, 1 )
            fprintf( fid, '%20f ', txtTable( l, 1 ) );   % Event time after correction (in sample counts)
            fprintf( fid, '%6d ',  txtTable( l, 2) );    % Antenna ID 
            iid = find(aid == txtTable( l, 2) );
            if size(iid,1)>0
                lev = c(iid(end),3);
            else
                lev = 0;
            end
            fprintf( fid, '%6d ',  lev+1 );    % Event number
            fprintf( fid, '%6d ',  cid+1);    % Coinc number
            fprintf( fid, '%8.2f ', txtTable( l, 5 ) );  % Delay with 1rst antenna (standard method)
            fprintf( fid, '%8.2f ', txtTable( l, 6 ) );  % Delay with 1rst antenna (intercorrelation method)
            fprintf( fid, '%6.2f ', txtTable( l, 7 ) );  % Average correlation for this antenna (intercorrelation method)
            fprintf( fid, '%6d ', txtTable( l, 8 ) );  %  Saturation flag
            fprintf( fid, '%6.2f ', txtTable( l, 9 ) );  % Amplitude
            fprintf( fid, '%6.3f ', 0 ); % Gain OBSOLETE
            fprintf( fid, '%6.3f ', 0);  % Amplitude method 1 (bline)
            fprintf( fid, '%6.3f ', 0);  % Amplitude method 2 (PSD)
            fprintf( fid, '\n' );
        end
        fclose(fid);
        clear txtTable;
        
        % DST
        % Setup
        if E=='1e17'
          RunSetup.Run=1017;
        elseif E=='5e17'
          RunSetup.Run=5017;
        elseif E=='8e16'
          RunSetup.Run=8016;
        elseif E=='5e16'
          RunSetup.Run=5016;
        elseif E=='2e17'
          RunSetup.Run=2017;
        end
        RunSetup.TotalEvt=sum([dstDet.Evt]);
        Struct.Setup = RunSetup;
        Struct.Setup.Det = dstDet;
        Struct.Setup.DetMalfunction=0;
        Struct.Setup.DetNoTrig=0;
        Struct.Setup.TotalCoinc = dstNbOfCoinc+1;
        Struct.Setup.RunTimeStart = 0;
        Struct.Setup.RunTimeStop = 0;
        Struct.Setup.InfosRun.TimeStart = 0;
        Struct.Setup.InfosRun.TimeStop = 0;
        
        % Coinc
        tag = dstCoinc.Det.Tag(end,:);
        dstCoinc.Det.TrigCor(end,tag==1) =  dstCoinc.Det.TrigCor(end,tag==1)-min(dstCoinc.Det.TrigCor(end,tag==1));
        dstCoinc.Det.TrigTime(end,tag==1) =  dstCoinc.Det.TrigTime(end,tag==1)-min(dstCoinc.Det.TrigTime(end,tag==1));        
        dstCoinc.Mult(end+1) = trig(1);
        dstCoinc.MultAnt(end+1) = trig(1);
        dstCoinc.IdCoinc(end+1) = cid+1;
        if dstCoinc.IdCoinc(end) ~= cid+1
            disp 'Error!'
            dstCoinc.IdCoinc
            pause
        end
        dstCoinc.IsShower(end+1) = 1;
        Struct.Coinc = dstCoinc;      
        
        % Simu
        dstSimu{end+1}.JobId = jobid;
        if jobid < 1000000
            dstSimu{end}.E = 5e17;
        elseif jobid>1000000
            dstSimu{end}.E = 8e16;
        end
        dstSimu{end}.Dir = [theta phi];
        dstSimu{end}.Core = [xCore yCore];
        dstSimu{end}.nAnt = ispres(1);
        dstSimu{end}.nTrig = trig(1);
        dstSimu{end}.Core = [xCore yCore];        
        Struct.Simu = dstSimu;       

        filename =sprintf(dst_filename,RunSetup.Run,1);
        save(filename,'Struct');
        display(sprintf('DST %s now saved to file.',filename));

     end

  %% Plots
  if DISPLAY >0
      maxi(find(maxi(:,1)==0),:) = [];
      maxi(:,2:7) = maxi(:,2:7)/1000;
      maxi
      %
      % Eta & Beta
      figure(11)
      subplot(2,1,1)
      plot(id,eta,'sk-','LineWidth',2,'MarkerFaceColor','k')
      hold on
      grid on
      plot(id,eta_v,'sk--','LineWidth',2)
      xlabel('Antenna Id', labelOpts{:})
      ylabel('Eta [deg]', labelOpts{:})
      ylim([0 90])
      subplot(2,1,2)
      plot(id,beta,'sk-','LineWidth',2,'MarkerFaceColor','k')
      hold on
      grid on
      plot(id,beta_v,'sk--','LineWidth',2)
      xlabel('Antenna Id', labelOpts{:})
      ylabel('Beta [deg]', labelOpts{:})
      ylim([0 90])
      %
      eta_m = mean(eta(eta>0));
      beta_m = mean(beta(beta>0));

      % Lateral distrib
      in = maxi(:,1);
      [a indant ] = intersect(idant,in);
      posDet = [x(indant) y(indant) z(indant)];
      posCore = (ones(length(in),1)*[xCore yCore zCore]);
      dX =  posDet - posCore;
      cp = cosd(phi); 
      sp = sind(phi);
      ct = cosd(theta); 
      st = sind(theta);
      u  = [ -sp*st, cp*st, ct ];  %% Warning!!!! This is different from usual since here x=WE & y=SN
      U  = ones( length(in), 1 )*u;
      dist = sqrt( sum( dX.^2, 2 ) - sum( dX.*U, 2 ).^2 ); %Distance from antennas to shower axis     

      norm = max(max(maxi(:,2:4)));
      for i = 1:size(tit,2)
          figure(i)
          set(i,'Name',tit{i},'NumberTitle','off')
          subplot(2,1,1)      
          plot(x,y,'k^','MarkerSize',8,'MarkerFaceColor','y')
          hold on
          grid on
          plot(xCore,yCore,'mh','MarkerSize',8,'MarkerFaceColor','m')
          xlabel('Easting [m]',labelOpts{:})
          ylabel('Northing [m]',labelOpts{:})
          %
          norm = max(max(maxi(:,i+1)));
          %norm = 1;
          for j=1:size(maxi,1)
            siz = max(1,10*log(10*maxi(j,i+1)/norm+1));
            plot(x(indant(j)),y(indant(j)),'o','MarkerSize',siz,'Linewidth',siz,'Color', col{i})
            text(x(indant(j))+20,y(indant(j)),sprintf('%d',maxi(j,1)))
          end
          plot(xCore,yCore,'mh','MarkerSize',8,'MarkerFaceColor','m')
          %
          xtraj = [-500 3000];
          ytraj = tand(phi_eva)*(xtraj-xCore)+yCore;
          plot(xtraj,ytraj,'--m')
          axis([min(xCore,min(x))-200 max(xCore,max(x))+200 min(yCore,min(y))-200 max(yCore,max(y))+200])
          %
          subplot(2,1,2)
          semilogy(dist,maxi(:,i+1),'sk','MarkerSize',8,'MarkerFaceColor','k')
          hold on
          for j=1:size(maxi,1)
            %ind = find(id == maxi(j,1));
            text(dist(j)+20,maxi(j,i+1),sprintf('%d',maxi(j,1)))
          end
          hold on
          grid on
          xlabel('Distance to shower axis [m]',labelOpts{:})
          if i<4
              ylabel('Amplitude [mV/m]',labelOpts{:})
          else
              ylabel('Voltage [mV]',labelOpts{:})
          end
          amp(i) = mean(maxi(dist<300,i+1));
      end  
      figure(1)
      figure(4)
      figure(11)
      pause
  end
  
  close all
  toc
  
  if DISPLAY
      figure(72)
      plot(x,y,'k^','MarkerSize',8,'MarkerFaceColor','y')
      hold on
      grid on
      %plot(coreAll(find(trig(:,1)>=5),1),coreAll(find(trig(:,1)>=5),2),'hg','MarkerFaceColor','g','MarkerSize',8);  % Triggered shower
      %plot(coreAll(find(trig(:,1)<5),1),coreAll(find(trig(:,1)<5),2),'hk','MarkerFaceColor','k')
      title(sprintf('Core positions for shower (%d,%d) degs',theta,phi))
      xlabel('Easting [m]',labelOpts{:})
      ylabel('Northing [m]',labelOpts{:})
  end
  
  
%%  
function [tvr AntDataRes] = computeDigitizedSignal(tt,AntData)
     SharedGlobals;
     %
     stepi = mean(diff(tt))*1e-9; %s % Time step for tt
     r = round(1/(stepi*FSAMPLING)); % Ratio between final time step and ini time step
     %
     %phase = r*rand(1,1); % Not needed (random anyway)
     ind = (1:floor(length(tt)/r))*r;
     AntDataRes = AntData(ind);
     tvr = tt(ind);
     %tvr = decimate(tt,r)+2.5;
     %tvr = (1:length(AntDataRes))*r*stepi*1e9; %ns

