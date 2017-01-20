function [senstot sensall dur] = ComputeAcceptance(nrun,sub,dstf)
% Compute efficiency skyplot for given run
% sensall: taking into account T6 (TagOut antennas for all directions)
% OMH 21/11/2012

SharedGlobals;
dur = 0;
senstot = zeros(1,360);
sensall = zeros(1,360);
%% Load dst
if ~exist('dstf')
    dstname = [DST_PATH sprintf(dst_filename,nrun,sub)];
    disp(sprintf('Loading dst %s',dstname))
    dst = load(dstname);
    Struct = dst.Struct;
else
    Struct = dstf;
end
SetupStruct = Struct.Setup;
Detectors = [SetupStruct.Det.Name];
CoincStruct = Struct.Coinc;
Mult = CoincStruct.MultAnt;
nCoincs = length(Mult);
display(sprintf('%d coincs in dst.',nCoincs));
if nCoincs == 0
  disp 'Skip.'
  return
end
time = max(CoincStruct.Det.UnixTime,[],2)/60;  % minutes
date_start = min(time);
time = time-date_start;
dur = ceil(max(time));
PhiP=Struct.Coinc.PlanRecons.Radio.Phi;
%PhiS=Struct.Coinc.SphRecons.Phi;
Tag=Struct.Coinc.Det.Tag;
Slope = Struct.Coinc.PlanRecons.Radio.SlopeDelay;
Chi2 = Struct.Coinc.PlanRecons.Radio.Chi2Delay;
Radius = Struct.Coinc.SphRecons.minDistSource;
TagTot = sum(Tag,1);
AntsIn = find(TagTot>100);  %Antenna alive if more than 100 trigs

%%
phi = cell(1,dur);
tag = cell(1,dur);
tagall = cell(1,dur);

T5_span = cutsettings.T5.time_span_T5_mn;

for i=1:dur
    if i/10==floor(i/10)
        disp(sprintf('Loop 1/2: %d / %d mns',i,dur))
    end
    %disp(sprintf('Now in slot %d-%d / %3.1f mns',i-1,i,max(time)));
    sel = find(time>i-1 & time<=i);
    if size(Radius,1)>1
        Radius = Radius';
    end
    selg =  intersect(sel,find(abs(Slope'-1)<0.1 & Chi2'<50 & Radius'>cutsettings.T5.DistCut));
    %disp(sprintf('%d coincs, %d distant & valid',length(sel),length(selg)))
    
    tagall{i} = Tag(sel,:);
    for k = i-T5_span:i+T5_span
        if k>0 & k<=dur
            phi{k} = [phi{k} PhiP(selg)]; 
            tag{k} = [tag{k} Tag(selg,:)'];
%              figure(1)
%              i
%              k
%              hist([phi{k}],100)
%              length([phi{k}])
%             pause
%             close(1);
        end
    end  
end

%%
senstot = zeros(1,360);
sensall = zeros(1,360);
for i=1:dur
    if i/10==floor(i/10)
        disp(sprintf('Loop 2/2: %d / %d mns',i,dur))
        ok = 1
    end

    %% Same antennas events
    ThisCoincs = [tagall{i}];
    ncoinc = min(150,size(ThisCoincs,1)); %Limit max number of coincs scanned to 200... otherwise too slow
    TagOut = zeros(1,size(ThisCoincs,2));
    disp(sprintf('%d/%d mns: %d coincs',i,dur,ncoinc))
    tic
    for j = 1:ncoinc-1
        for k = j+1:ncoinc
          %com = intersect(find(ThisCoincs(j,:)>0),find(ThisCoincs(k,:)>0));  %intersect is too slow
          com = find(ThisCoincs(j,:)+ThisCoincs(k,:)==2);
%           ncoinc
%           j
%           k
%          length(com)
          if length(com)>=cutsettings.T6.NbComAntT6
            if HARD  % HARD cuts: only  2 events in common are enough
               TagOut(com) = 1;  % Tag all antennas common to these 2 events
            else  % SOFT cuts: need 3 events in common 
                for l = k+1:ncoinc % Look for 3rd event
                    [com3 ia] = find(ThisCoincs(j,com)+ThisCoincs(l,com)==2);
                    if length(com3)>=cutsettings.T6.NbComAntT6  % 3rd event was found.
                        TagOut(com(ia)) = 1;  % Tag all antennas common to these 3 events
                    end
                end
            end
          end
          %pause
        end
    end
    toc
    ok = 2
    %pause
    %
    sens = zeros(1,360);
    senst = zeros(1,360);   
    phik = phi{i};
    tagk = tag{i}';
%     figure(1)
%     hist([phi{i}],100)
%     i
%     length([phi{i}])
%     pause
%     close(1)
    %% Same direction events
    for phii=1:360
        dsens = zeros(1,360);
        dsens2 = zeros(1,360);
        
        phimin = phii-cutsettings.T5.PhiCut;
        phimax = phii+cutsettings.T5.PhiCut; 
        if phimin<0 
            sel = find( (phik<phimax & phik>0) | (phik>360+phimin & phik<360) );
        elseif phimax>360
            sel = find( (phik>phimin &  phik<=360) | (phik>0 & phik<phimax-360) );      
        else
            sel = find(phik<phimax & phik>=phimin);
        end
        ntrig = sum(tagk(sel,:),1);
%        [Detectors' ntrig' TagOut']        
        nin = sum(ntrig>0);
        r = 1-nin/length(AntsIn);    
        nin2 = min(length(AntsIn),sum((ntrig+TagOut)>0));
        r2 = 1-nin2/length(AntsIn);
        if phimin<=0 
            dsens(1:phimax) = r;
            dsens(360+phimin:360) = r;
            dsens2(1:phimax) = r2;
            dsens2(360+phimin:360) = r2;
        elseif phimax>360
            dsens(1:phimax-360) = r;
            dsens(phimin:360) = r;
            dsens2(1:phimax-360) = r2;
            dsens2(phimin:360) = r2;           
        else
            dsens(phimin:phimax) = r;
            dsens2(phimin:phimax) = r2;            
        end
%         if nin2>0
%         disp(sprintf('[%d-%d deg]: %d antennas triggered in slice %d-%d /%3.1f mins',phimin,phimax,nin,i-1,i,dur))        
%          nin2
%          nantsin=length(AntsIn)
%          nout=sum((ntrig+TagOut)>0)
%          r2
%          dsens2
%          pause
%         end
        sens = sens + dsens;
        senst = senst + dsens2;
%         if nin>0
%         figure;
%         plot(sens,'k-')
%         hold on
%         plot(senst,'r-')
%         pause
%         end
%        close all
    end
    norm = 2*cutsettings.T5.PhiCut+1;
    sens = sens/norm;
    senst = senst/norm;
    
%         figure;
%         plot(sens,'k-')
%         hold on
%         plot(senst,'r-')
%         pause
%         close all
    
    senstot = senstot + sens;
    sensall = sensall + senst;
end
senstot = senstot./dur;
sensall = sensall./dur;

%% Save to file
filename = [ACC_PATH sprintf('Acceptance_R%d_%d.mat',nrun,sub)];
save(filename,'senstot','sensall','dur')

clear dst;

%% Plot
% figure(10);
% set(10,'Name',sprintf('Acceptance R%d',nrun),'NumberTitle','off')
% plot(1:360,senstot,'k','LineWidth',2)
% hold on
% plot(1:360,sensall,'r','LineWidth',2)
% xlabel('Azimuth angle [deg]',labelOpts{:})
% ylabel('Relative acceptance',labelOpts{:})
% xlim([0 360])
% %ylim([0 1])
% grid on
% hold on

