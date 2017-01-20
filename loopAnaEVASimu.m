function loopAnaEVASimu(E)
% loop on AnaEVASimu
% OMH 24/03/2014

SharedGlobals;

if ~exist('E')
  E = '1e17';
end
disp(sprintf('\nEnergy = %s eV.', E))

%necfolder = '../data/simu/trend-50/1*';
necfolder = [SIMU_PATH E '/voltages/*']

%fl = ls(necfolder)
bb = dir(necfolder);
for i=1:size(bb,1)
  fl(i,1:length(bb(i).name))=bb(i).name;
end

jobids=[];
for k = 1:size(fl,1)
   c =  strread(fl(k,:),'%s','delimiter','_'); 
   % Cut folder name into pieaces separeted by '_' as stdard name is jobID_Xcore_Ycore.
   if size(c,1)>1  & size(c{1},2)>6 % Standard name
      jobids(end+1) = str2num(c{1});
   end
   %if size(str2num(fl(k,1:sjobid)),1)>0
   %    jobids(end+1)=str2num(fl(k,1:sjobid));
       %Workaround
   % end
end
jobids = unique(jobids);  % remove redondant jobs
clear fl;
trigAll = [];
presAll = [];
polAll = [];
coreAll =  [];

thetaAll = zeros(length(jobids),1);
phiAll = zeros(length(jobids),1);
res = zeros(length(jobids),3);

%% Load skymap
filename = sprintf('simuSkymap_%s.mat',E);
if fopen(filename)>0
    a = load(filename);
    done = a.coreAll(:,[1 4 5]);
else
  done = [];
end

%jobids = 8630099;

%% Loop on showers
%for k = 1:length(jobids)
for k = 1:550
%length(jobids)
    if jobids(k)==8665868 | jobids(k)==8629961 
        continue
    end

    necfolder = sprintf('%d*',jobids(k));
    %necfolder = ['../data/simu/trend-50/' necfolder];
    necfolder = [SIMU_PATH E '/voltages/' necfolder];
    %fl = ls(necfolder)
    %clear fl;
    bb = dir(necfolder);
    if size(bb,1)<1
      disp(sprintf('No folders %s',necfolder))
      bb
      clear necfolder;
      clear bb
      %pause
      %continue
      break
    end
    for i=1:size(bb,1)
       fl(i,1:length(bb(i).name))=bb(i).name;     
    end
    fl(size(bb,1)+1:end,:)=[];

    disp(sprintf('*** Shower %d/%d \nNow processing %d folders associated to shower %d:',k,length(jobids),size(fl,1),jobids(k)))
    
    %% Loop on core positions
    thisTrig = zeros(1,3);
    thisPres = zeros(1,3);
    for i = 1:size(fl,1)
      %fl(i,:)
      %size(fl)
      %
      truename= strread(fl(i,:),'%s','delimiter',' ');
      %folder = sprintf('../data/simu/trend-50/%s/',truename{1});
      folder = sprintf('%s/%s/voltages/%s/',SIMU_PATH,E,truename{1});
      [ast theta phi_eva xCores zst] = strread(folder,'%s%d%d%s%s','delimiter','_');
      phi = mod(phi_eva-90,360);
      %
      xCore = str2num(xCores{:});
      if size(xCore,1) == 0
          disp(sprintf('*** Core position %d/%d: Faulty xCore value. Skip this core position.',i,size(fl,1)))
          continue
      end
      %
      zzst = zst{1,:};
      yCore = str2num(zzst(1:end-1)); 
      if size(yCore,1) == 0
          disp(sprintf('*** Core position %d/%d: Faulty yCore value. Skip this core position.',i,size(fl,1)))
          continue
      end
      zCore = REFALT;

      %
      % Check if this core pos was processed already...
      ind = find(sum(ismember(done,[jobids(k) xCore yCore]),2)==3);
      if size(ind,1)>0 % Already processed
          disp(sprintf('Core position (%dm,%dm) already present at line %d in %s. Skip treatment.',xCore,yCore,ind,filename))
          thisTrig = thisTrig+(a.trigAll(ind,:)>=5);
          thisPres = thisPres+(a.presAll(ind,:)>=5);       
          theta = a.coreAll(ind,2);
          phi = a.coreAll(ind,3);
          trigAll(end+1,:) = a.trigAll(ind,:);
          presAll(end+1,:) = a.presAll(ind,:);
          polAll(end+1,:) = a.polAll(ind,:);
          coreAll(end+1,:) = a.coreAll(ind,:);   
          continue
      end

      if theta>80 
          disp(sprintf('Inclined shower. Skip.'))
          continue
      end
      if theta<60 
          disp(sprintf('Already processed in previous skymap.mat Skip.'))
          continue
      end
      %if phi>90 
      %    disp(sprintf('Shower outside Junhua s sector. Skip.'))
      %    continue
      %end
      
      %% Check trajectory
      if checkTraj(theta,phi,xCore,yCore)==0
          disp(sprintf('Core position %d (%d,%d)m :shower trajectory below ground level!\nImpose detection = 0.',i,xCore,yCore))
          bb=dir([folder '/*ew.txt']);
          nfew = size(bb,1);
          bb = dir([folder '*ns.txt']);
          nfns = size(bb,1);
          bb = dir([folder '*vert.txt']);
          nfv = size(bb,1);
          %
          trig = [0 0 0];
          pol = zeros(1,7);
          ispres = [size(nfew,1) size(nfns,1) size(nfv,1)];
      else
          disp(sprintf('*** Core position %d/%d: now calling AnaEVASimu(%d,%s,%d,%d)...',i,size(fl,1),jobids(k),E,xCore,yCore))
          [trig ispres theta phi pol] = AnaEVASimu(jobids(k),E,xCore,yCore);
      end      
      
      trigAll(end+1,:) = trig;  % Number of antennas triggered / channel
      presAll(end+1,:) = ispres;  % Number of antennas present / channel
      polAll(end+1,:) = pol;  % Number of antennas present / channel
      coreAll(end+1,:) = [jobids(k) theta phi xCore yCore];
      thisTrig = thisTrig+(trig>=5);
      thisPres = thisPres+(ispres>=5);
    end  % Loop on Core pos
    clear fl;
    clear bb;
    fclose all;
   
    thetaAll(k) = theta;
    phiAll(k) = phi;
    res(k,:) = thisTrig./thisPres;
    disp ' '
    disp(sprintf('Fraction of showers detected for (%d, %d)deg: %3.1f pc for EW polar, %3.1f pc for NS polar.',theta, phi,res(k,1)*100,res(k,2)*100))    
    %
    allRes = [jobids' thetaAll phiAll res];
    save(sprintf('simuSkymap_%s.mat',E),'allRes','trigAll','presAll','coreAll','polAll');
end
allRes = [jobids' thetaAll phiAll res];
allRes;
trigAll;
save(sprintf('simuSkymap_%s.mat',E),'allRes','trigAll','presAll','coreAll','polAll');
