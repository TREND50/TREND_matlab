function [antout lt me mde] = CheckSimuSig(jobid,E,xCore,yCore)
% Checks if Efield signals are OK (window duration)
% Taken from AnaEVASimu
% OMH 04/02/2016

  SharedGlobals;
  CC = 1;
  DISPLAY = 0;
  
  col = {'b','r','g','b','r','g'};
  tit = {'E_{EW}','E_{NS}', 'E_{vert}','EW antenna voltage','NS antenna voltage','Vert antenna voltage'};
  f = {'efield.txt','ew.txt','ns.txt','vert.txt'};
  
  %% Set path & load antenna coordinates    
  FREQMIN = 50;  % 50-100MHz bandwidth for TREND
  if CC == 0
    SPATH = '../data/simu/TREND50/';
    RPATH = SPATH;  % Result path
    coord = load([RPATH 'coord_antennas_TREND50.txt']);
  else
    SPATH = '/sps/hep/trend/trend-50/';
    RPATH = './TREND50/';  % Result path
   coord = load([RPATH 'coord_antennas_TREND50.txt']);
  end
  
  idant = coord(:,1);
  id = idant;
  x = coord(:,2);
  y = coord(:,3);
  z = coord(:,4);
  
  amp = zeros(1,6);
  pol = zeros(1,7);
  trig = zeros(1,3);
  ispres = zeros(1,3);
  theta = 0;
  phi = 0;
  antout = [];
  lt = [];
  me = [];
  mde = [];
  
  necfolder = sprintf('/%d*%d_%d*',jobid,xCore,yCore);
  necfolder = [SPATH E '/voltages/' necfolder];  
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
  %
  truename= strread(fl(:),'%s','delimiter',' ');
  clear fl;
  clear dd;
  fclose all;
  folder = sprintf('%s/%s/voltages/%s/',SPATH,E,truename{1});
  [ast theta phi_eva xCore zst] = strread(folder,'%s%d%d%d%s','delimiter','_');
  phi = mod(phi_eva-90,360);
  zCore = REFALT;

  tic
  disp(sprintf('Zen = %d deg, Az = %d deg, core position = (%dm, %dm)',theta,phi,xCore,yCore))
  maxi = zeros(length(idant),10);
  txtTable = zeros(1,12);
  %ok = 0;

  %% Loop on antennas  
  for i = 1:length(idant)
      ind = find(idant==id(i));
      %disp(sprintf('Antenna %d: (%3.0f,%3.0f,%3.0f)',idant(ind),x(ind),y(ind),z(ind)))

      filename = sprintf('a%d_%s',idant(i),f{1});
      filename = [folder filename];
      if fopen(filename)<0
          %disp(sprintf('No file %s. Skip.',filename))
          continue
      end
        
      %% Loop on channels
      for j=1:3
        % Loop on fiel=ds, with:
        % 1: Efield
        % 2: EW
        % 3: NS
        % 4: Vert

        %% Load files
        filename = sprintf('a%d_%s',idant(i),f{j});
        filename = [folder filename];
        if fopen(filename)<0
          %disp(sprintf('No file %s. Skip.',filename))
          break
        end
        m = load(filename);  
        if size(m,1) == 0
          %disp(sprintf('File %s is empty. Skip.',filename))
          continue
        end
        %           
        if j == 1  % Efield
          te = m(:,1)*1e-6; % s
          xe = m(:,2); % Ex component [muV/m] (EW)
          ye = m(:,3); % Ey component [muV/m] (NS)
          ze = m(:,4); % Ez component [muV/m] (Vert)
          %
          dxe = [0; diff(xe)];
          dye = [0; diff(ye)];
          dze = [0; diff(ye)];
          te = te-te(1);
          tens = te*1e9;  %ns
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
          dv(j-1,:) = [0; diff(vs)];
          t = t-t(1);
          tns = t*1e9;  %ns
        end
        
        %% Plot
        if DISPLAY
              if j == 1
                  figure(1)
                  subplot(2,1,1)
                  plot(tens,xe,'b','LineWidth',2)
                  hold on
                  grid on
                  plot(tens,ye,'r','LineWidth',2)
                  plot(tens,ze,'g','LineWidth',2)                  
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Efield [muV/m]', labelOpts{:})
                  title(sprintf('Antenna %d - Efield',idant(i)),labelOpts{:})
                  hold off
                  subplot(2,1,2)
                  plot(tens,dxe,'b','LineWidth',2)
                  hold on
                  grid on
                  plot(tens,dye,'r','LineWidth',2)
                  plot(tens,dze,'g','LineWidth',2)                  
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Diff Efield [muV/m/ns]', labelOpts{:})
                  hold off
              else
                  figure(2)
                  subplot(2,1,1)
                  plot(tns(tns<200),v(j-1,tns<200),'Color',col{j-1},'LineWidth',2)
                  hold on
                  grid on
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Voltage [muV]', labelOpts{:})
                  title(sprintf('Antenna %d - %s',idant(i),tit{j+2}),labelOpts{:})
                  subplot(2,1,2)
                  plot(tns(tns<200),dv(j-1,tns<200),'Color',col{j-1},'LineWidth',2)
                  hold on
                  grid on
                  xlabel('t [ns]', labelOpts{:})
                  ylabel('Diff Voltage [muV/ns]', labelOpts{:})
                  title(sprintf('Antenna %d - %s',idant(i),tit{j+2}),labelOpts{:})
              end
        end
      end   % loop on channels         
      antout(end+1) = idant(i);
      lt(end+1) = tens(end);
      me(end+1) = abs(xe(end)./max(abs(xe))); 
      mde(end+1) = abs(dxe(end)./max(abs(dxe))); 
      %tns(end) 
      %max(dv(:,end))
  
      if DISPLAY
        figure(2)
        subplot(2,1,1)
        hold off
        subplot(2,1,2)
        hold off
        pause
      end
      
      fclose all;
      clear v;
      clear dv;
      clear v_r;
      clear v_th;
      clear v_ph;

  end  % loop on antennas
