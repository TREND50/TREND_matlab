function [EventTimeTable,Struct]=EmptySignalFilter(EventTimeTable,Struct)

SharedGlobals;

ant=[Struct.Setup.Det.Name];
nrun=Struct.Setup.Run;

EmptyFlag=zeros(size(EventTimeTable,1),1);

for i=1:length(ant)
    
    % Events for the detector
    ind_event=find(EventTimeTable(:,2)==ant(i));
    % Open file
    fd = OpenFileData(nrun,ant(i));
    
    if fd==-1
        display(sprintf('No data file for antenna %d',ant(i)))
        continue
    else
        for j=1:length(ind_event)
            
            % Find concerned event in raw data
            nevt=EventTimeTable(ind_event(j),3);
            fseek(fd,ibuff*(nevt-1),'bof');
            DataEvt=double(fread(fd,ibuff,'uint8'));
            
            if max(DataEvt)<129 & min(DataEvt)>109 % Check if signal not too weak (+/- 10 around mu)
                S=[];
            else           
                % Sigma and Mu calculation
                KLOW=[1:(round(length(DataEvt)/2)-100)];
                indmax=find(max(DataEvt(KLOW))==DataEvt(KLOW));
                if DataEvt(KLOW(indmax))>155 & indmax>150
                    KLOW=[1:indmax-10];
                end;
            
                if EmptySignalTotalFilt==1
                    [ y, x ] = hist( DataEvt( KLOW ), [ 0:255 ] );
                    [ ym, km ] = max( y );
                    xm  = x( km );
                    Dx  = unique( abs( x( y > 0.5*ym ) - xm ) ); % Select samples within Width at Half Full Height (WHFH).
                    
                    if length(Dx)==0 |  Dx==0
                        mu = mean(DataEvt(KLOW));
	                sig=std(DataEvt(KLOW));
                    else
                        n = round( Dx(end)*2.5577 ); % 3 sigmas width, corrected from WHFH and rounded.  
                        K   = km + [ -n:1:n ];
                        K   = K( K <= 256 & K >= 1 );
                        xk  = x( K );
                        pk  = y( K )/sum( y( K ) );
                        mu  = sum( xk.*pk );                 % Noise mean value
                        sig = sqrt( sum( ( xk - mu ).^2.*pk ) ); % Noise sigma
                    end
%                     threshold*sig
%                     length(find(abs(DataEvt-mu)>=threshold*sig))
%                     figure,plot(DataEvt)
%                     pause;
%                     close all
                else
                    mu=mean( DataEvt(KLOW));
                    sig=std( DataEvt(KLOW));
                end;
                     
                % Look for signal upper than threshold
                Win=(abs(DataEvt-mu)>=threshold*sig)';
                %Win=Agglomerate(Win,granularity);
                %Win =Agglomerate(Win,granularity,1);
                S =find(Win);
            end;
            
            if isempty(S) | isempty(find(S>400 & S<600))
                EmptyFlag(ind_event(j))=1;
            end;
            
            clear DataEvt KLOW mu sig Win S
            
        end
        
    end

display(sprintf('Treatment done for antenna %d',ant(i)))

end

% Remove empty signal from EventTimeTable and save informations
indok=find(EmptyFlag==0);
indempty=find(EmptyFlag==1);

EmptySignalInfo=EventTimeTable(indempty,:);
Struct.Setup.InfosRun.EmptySignals=EmptySignalInfo;

EventTimeTable=EventTimeTable(indok,:);
display(sprintf('%d empty events has been erased (%d percent)',length(indempty),100*length(indempty)/(length(indempty)+length(indok))))
                
            
            
