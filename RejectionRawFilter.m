function [TagFlag,FiltResult]=RejectionRawFilter(Carac,settings)
% Returns 0 if signal rejected, 1 if OK

SharedGlobals;

% Initial value of tagged antenna
TagFlag=1;

% Cut parameters
Threshold = settings.threshold;
Max_out_ToT   = settings.max_out_ToT;
Max_block_ToT   = settings.max_block_ToT;
Max_pretrig_ToT = settings.max_pretrig_ToT;

% Boxes characteristics
nbox=cell2mat(Carac(2));
box_start=cell2mat(Carac(4));
box_end=cell2mat(Carac(5));
box_length=cell2mat(Carac(6));
box_amp=cell2mat(Carac(7));
box_dt=cell2mat(Carac(8));

FiltResult=zeros(7,2);
% 1 : ToT 
% 2 : ToT before trigger (100ns)
% 3 : Central Box?  OUT
% 4 : ToT Central Box duration;

% nbox
% box_start
% box_amp
% Threshold

% Get central box
box_mid = (box_start+box_end)/2;
[a indcenter]=min(abs(box_mid-ibuff/2)); 
% Get other boxes
indout = 1:nbox;
indout(indcenter) = []; 

%% Test5: Large pulse within 1mus before central box
if nbox==0
    FiltResult(5,1)=nbox;
    FiltResult(5,2)=1;
    TagFlag=0;
    %disp 'No box'
    %pause
    return
end


%% Test2: Time before center
TotalPreTrig = (ibuff/2-box_start(indcenter))/FSAMPLING;
if TotalPreTrig > Max_pretrig_ToT
     FiltResult(2,1)=TotalPreTrig;
     FiltResult(2,2)=1;
     TagFlag=0;
%     disp 'Time before trigger out'
%      TotalPreTrig
%      Max_pretrig_ToT
%      pause
     return
end;

%% Test3: Central box duration
if box_length(indcenter) > Max_block_ToT
    FiltResult(3,1)=box_length(indcenter);
    FiltResult(3,2)=1;
    TagFlag=0;
%     disp 'Central box duration out'
%     box_length(indcenter)
%     Max_block_ToT
%     pause
    return
end;

%% Test1: ToT outside central box
if nbox>1
    time_over_threshold = sum(box_length(indout));
    if time_over_threshold > Max_out_ToT
        FiltResult(1,1)=time_over_threshold;
        FiltResult(1,2)=1;
        TagFlag=0;
%         disp 'Tot oustide box out'
%         time_over_threshold
%         Max_out_ToT
        return
    end
end

%% Test4: Large pulse within 1mus before central box
if nbox>1
    big = find(box_amp(indout)> Threshold);
    close = find(box_start(indout)<ibuff/2 & abs(box_mid(indout)-ibuff/2)<1e-6*FSAMPLING);
    large = find(box_length(indout)>50e-9);
    closebig = intersect(big,close);
    bad = intersect(closebig,large);
    %close
    %big
    if length(bad) > max_out
        FiltResult(4,1)=length(bad);
        FiltResult(4,2)=1;
        TagFlag=0;
        %disp 'Strong signal within 1mus'
        %pause
        return
    end
end


%% OBSOLETE
% % Test4: "Bounce" ratio (box with max amp very small in comparison with the
% % central box)
% 
% K1=(box_amp < BOUNCE);
% 
% 
% 
% 
% 
% % Cut 1: boxes multiplicity and total ToT
% 
% if nbox > Max_multiplicty 
%     nboxtemp=0;
%     FiltResult(1,2)=1;
% elseif time_over_threshold > Max_total_ToT
%     nboxtemp=0;
%     FiltResult(2,2)=1;
% elseif (nbox>0)
%     
%     % Cut 2 : bounce ratio, "too long" blocks, and "too close" blocks
%     ind=find(box_start<=(ibuff/2));
%     indcenter=ind(find( abs(box_start(ind)-(ibuff/2))==min(abs(box_start(ind)-(ibuff/2)))));
%     if isempty(indcenter)
%         indcenter=1;
%     end;
%     BOUNCE=Bounce_ratio*box_amp(indcenter);
%     K1=(box_amp < BOUNCE);
%     FiltResult(3,1)=sum(K1);
%     K2=(box_dt < Inhibit_window);
%     FiltResult(4,1)=sum(K2);
%     K3=(box_length > Max_block_ToT);
%     FiltResult(5,1)=sum(K3);
%     
%     K = ( box_amp >= BOUNCE ) & ...
%     ( box_length <= Max_block_ToT ); % & ...
%     %( box_dt >= Inhibit_window );
%     
%     nboxtemp = sum( K );
%     if nboxtemp==0
%         FiltResult(3,2)=1;
%     end
%     boxstarttemp=box_start( K );
%     boxendtemp=box_end( K );
% else
%     nboxtemp=0;
%     FiltResult(6,2)=1;
% end;
%      
% % box containing the trigger signal
% if TrigTimeBox
%     if nboxtemp>0
%         indtrig=find(boxstarttemp<=(ibuff/2) & boxendtemp>(ibuff/2));
% %         boxstarttemp
% %         boxendtemp
%         if isempty(indtrig)
%             TagFlag=0;
%             FiltResult(7,2)=1;
%         end;
%     else
%         TagFlag=0;      
%     end;
% else
%     if nboxtemp==0
%         TagFlag=0;
%     end;
% end;
