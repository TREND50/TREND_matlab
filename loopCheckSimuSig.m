function loopCheckSimuSig(E)
% Taken from loopAnaEVASimu
% OMH 04/02/2016

SharedGlobals;
CC = 1;

    if CC == 0
      SPATH = '../data/simu/TREND50/';
      RPATH = SPATH;  % Result path
    else
      SPATH = '/sps/hep/trend/trend-50';
      RPATH = './TREND50/';  % Result path
    end

if ~exist('E')
  E = '1e17';
end
disp(sprintf('\nEnergy = %s eV.', E))

necfolder = [SIMU_PATH E '/voltages/*']
%necfolder = [SPATH sprintf('/%s',E) '/voltages/*']

bb = dir(necfolder);
%fl = [];
for i=1:size(bb,1)
  fl(i,1:length(bb(i).name))=bb(i).name;
end

jobids=[];  
% Ugly trick to read the right number of digits in the name...
if E == '5e17'
  sjobid = 7;
else
  sjobid = 8;
end

resAll = [];
for k = 1:size(fl,1)
    c =  strread(fl(k,:),'%s','delimiter','_'); 
    % Cut folder name into pieaces separeted by '_' as stdard name is jobID_Xcore_Ycore.
    if size(c,1)>1  & size(c{1},2)>5  % Standard name
       jobids(end+1) = str2num(c{1});
    end
    if size(str2num(fl(k,1:sjobid)),1)>0
        jobids(end+1)=str2num(fl(k,1:sjobid));
    end
end
jobids = unique(jobids);  % remove redondant jobs
clear fl;

  
%% Loop on showers
for k = 1:length(jobids)
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
    for i = 1:size(fl,1)
      %
      truename= strread(fl(i,:),'%s','delimiter',' ');
      
      %folder = sprintf('../data/simu/trend-50/%s/',truename{1});
      folder = sprintf('%s/%s/%s/',SIMU_PATH,E,truename{1});
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
%       % Check if this core pos was processed already...
%       ind = find(sum(ismember(done,[jobids(k) xCore yCore]),2)==3);
%       if size(ind,1)>0 % Already processed
%           disp(sprintf('Core position (%dm,%dm) already present at line %d in %s. Skip treatment.',xCore,yCore,ind,filename))
%           thisTrig = thisTrig+(a.trigAll(ind,:)>=5);
%           thisPres = thisPres+(a.presAll(ind,:)>=5);       
%           theta = a.coreAll(ind,2);
%           phi = a.coreAll(ind,3);
%           trigAll(end+1,:) = a.trigAll(ind,:);
%           presAll(end+1,:) = a.presAll(ind,:);
%           polAll(end+1,:) = a.polAll(ind,:);
%           coreAll(end+1,:) = a.coreAll(ind,:);   
%           continue
%       end

      if theta>80 
          disp(sprintf('Inclined shower. Skip.'))
          continue
      end
      
      disp(sprintf('*** Core position %d/%d: now calling CheckSimuSig(%d,%s,%d,%d)...',i,size(fl,1),jobids(k),E,xCore,yCore))
      [ant lt me mde] = CheckSimuSig(jobids(k),E,xCore,yCore);
      posAnt = GetAntennaPosition(ant,0);
      if size(posAnt,2)>0
         [antid dist] = computeShowerDistance(theta,phi,[xCore,yCore],[ant' posAnt]);
         on = ones(length(antid),1);
         resAll(end+1:end+length(antid),:) = [jobids(k)*on theta*on phi*on xCore*on yCore*on antid dist lt' abs(me)' abs(mde)'];
         save('checkEVASignals.mat','resAll');
       else
          ant
       end
      %pause
    end  % Loop on Core pos
    clear fl;
    clear bb;
    fclose all; 
end

save('checkEVASignals.mat','resAll');
