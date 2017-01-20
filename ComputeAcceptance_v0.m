function [senstot sensall dur] = ComputeAcceptance(nrun,sub,dstf)
% Compute efficiency skyplot for given run
% OMH 21/11/2012

SharedGlobals;

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
        
%%
phi = cell(1,dur);
tag = cell(1,dur);
tagall = cell(1,dur);

for i=1:dur
    if i/10==floor(i/10)
        disp(sprintf('Loop 1/2: %d / %d mns',i,dur))
    end
    %disp(sprintf('Now in slot %d-%d / %3.1f mns',i-1,i,max(time)));
    sel = find(time>i-1 & time<=i);
    selg =  intersect(sel,find(abs(Slope'-1)<0.1 & Chi2'<50 & Radius'>500));
    %disp(sprintf('%d coincs, %d distant & valid',length(sel),length(selg)))
    
    tagall{i} = Tag(sel,:);
    for k = i-5:i+5
        if k>0 & k<=dur
            phi{k} = [phi{k} PhiP(selg)]; 
            tag{k} = [tag{k} Tag(selg,:)'];
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
          %com =
          %intersect(find(ThisCoincs(j,:)>0),find(ThisCoincs(k,:)>0));  %intersect is too slow
          com = find(ThisCoincs(j,:)+ThisCoincs(k,:)==2);
          if length(com)>=4  % Look for 3rd event
            for l = k+1:ncoinc
                [com3 ia] = find(ThisCoincs(j,com)+ThisCoincs(l,com)==2);
                if length(com3)>=4  % 3rd event was found.
                    TagOut(com(ia)) = 1;  % Tag all antennas common to these 3 events
                end
            end
          end
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
    %% Same direction events
    for phii=1:360
        dsens = zeros(1,360);
        dsens2 = zeros(1,360);
        
        phimin = phii-10;
        phimax = phii+10; 
        if phimin<0 
            sel = find( (phik<phimax & phik>0) | (phik>360+phimin & phik<360) );
        elseif phimax>360
            sel = find( (phik>phimin &  phik<=360) | (phik>0 & phik<phimax-360) );      
        else
            sel = find(phik<phimax & phik>=phimin);
        end
        ntrig = sum(tagk(sel,:),1);
        [Detectors' ntrig' TagOut'];        
        nin = sum(ntrig>0);
        r = 1-nin/50;    
        nin2 = sum((ntrig+TagOut)>0);
        r2 = 1-nin2/50;
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
        
        %disp(sprintf('[%d-%d deg]: %d antennas triggered in slice %d-%d /%3.1f mins',phimin,phimax,nin,i-1,i,dur))        
        sens = sens + dsens;
        senst = senst + dsens2;
%         figure;
%         plot(sens,'k-')
%         hold on
%         plot(senst,'r-')
%         pause
%         close all
    end
    sens = sens/21;
    senst = senst/21;
    
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
filename = sprintf('Acceptance_R%d_%d.mat',nrun,sub);
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

