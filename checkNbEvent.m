function totev = checkNbEvent(nrun)
% Counts nb of events per antenna for 2 different DST sets
% OMH 09/03/2014

SharedGlobals;

totev = zeros(2,length(ALLDETS));
%dstname = [DST_PATH sprintf('dst%d_%d.mat',nrun,1)];
dir = {'dst102012','dst102014'};

for i = 1:2
  dstname = ['/sps/hep/trend/' sprintf('%s/dst%d_%d.mat',dir{i},nrun,1)];
  if fopen(dstname)<0
        disp 'No dst available. Abort'
        return
  end
 disp 'Loading dst...'
  d = load(dstname);
  disp 'Done.'
  id=[d.Struct.Setup.Det.Name];
  [c,ia, ib] = intersect(id,ALLDETS);
  nev = [d.Struct.Setup.Det.Evt];
  totev(i,ib) = nev;
end
  
totev'

%plot(ALLDETS,totev(1,:)-totev(2,:),'+')
