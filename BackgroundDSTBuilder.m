function []=BackgroundDSTBuilder(nrun)

SharedGlobals;

% Background parameters
nbevent=1024;
size_event=1024;
time_buff=4;
FFTIN = 0;

if CC==1
  unix(sprintf('/sps/trend/fetchBack.sh %d',nrun))
end

[nbloop,ant,big]=ExtractInfosBckgr(nrun);
disp(sprintf('BackgroundDSTBuilder: performing background data analysis for R%d (%d antennas and %d loops)',nrun,length(ant),max(nbloop)))
Sigma=zeros(max(nbloop),length(ant));
Mu=zeros(max(nbloop),length(ant));
Time=zeros(max(nbloop),length(ant));
Spectrum=cell(max(nbloop),length(ant));

%big = 1;
% Loop on antennas
for i=1:length(ant)
        disp(sprintf('Antenna %d',ant(i)))
        %% Extraction du run
        fd = OpenFileBackground(nrun,ant(i));
        ft = OpenFileBackgroundTime(nrun,ant(i));
        if big==0
            data=fread(fd);
        end;
        time=fread(ft,inf,'uint32');            
        
        %% Loop on background loops
        for j=1:nbloop(i)
            if round(j/10)==j/10
                disp(sprintf('BackgroundDSTBuilder: analysing loop %d...',j))
            end
            % Event extraction
            start_evt=(j-1)*(nbevent*size_event)+1;
            if big==0
                RawEvt=data(start_evt:start_evt+nbevent*size_event-1);
            else
                fseek(fd,start_evt,'bof');
                RawEvt=fread(fd,nbevent*size_event);
            end;
            
            if ~isempty(RawEvt)                
                %% Signal filtering
                time_scale=[0:5e-9:(nbevent*size_event)*5e-9];
                Evt = mean(RawEvt)+PassBand(RawEvt,time_scale,FREQMIN,FREQMAX);
                Time(j,i)=time((j-1)*time_buff+1);

                %% Sigma/mu calculation (Valentin)
                [ y, x ] = hist( Evt, [ 0:255 ] );
                [ ym, km ] = max( y );
                xm  = x( km );
                Dx  = unique( abs( x( y > 0.5*ym ) - xm ) ); % Select samples within Width at Half Full Height (WHFH).
                n = round( Dx(end)*2.5577 ); % 3 sigmas width, corrected from WHFH and rounded.

                K   = km + [ -n:1:n ];
                K   = K( K <= 256 & K >= 1 ); 
                xk  = x( K );
                pk  = y( K )/sum( y( K ) ); 
                Mu(j,i)  = sum( xk.*pk );                     % Noise mean value
                Sigma(j,i) = sqrt( sum( ( xk - Mu(j,i) ).^2.*pk ) ); % Noise sigma

                if FFTIN
                %% Frequency spectrum (very slow)
                sumM=zeros(size_event,1);
                fleche_temps=[0:5e-9:(size_event-1)*5e-9];
                for k=1:nbevent
                    if size(Evt)<1024*k
                        disp(sprintf('Problem with event size for event %d... Skipping it.',k))
                        continue
                    end
                    DataEvt=Evt((k-1)*size_event+1:k*size_event);
                    [f,X,M]=tfou(DataEvt,fleche_temps);
                    sumM=sumM+M;
                end;             
                ind_freq=find(f>=FREQMIN);
                ind_freq=ind_freq(1:end-1);
                freq=f(ind_freq);
                Spectrum(j,i)={M(ind_freq)/nbevent};
                if DISPLAY
                    figure,plot(time_scale,Evt)
                    figure,plot(M(ind_freq)/nbevent)
                    figure
                    h = bar( x, y, 1.0 ); hold on;
                    set( h, 'FaceColor', 'none', 'EdgeColor', 'k' );
                    plot( x, ym*exp( -0.5*( x - mu ).^2/sig^2 ), 'r', 'LineWidth', 2 );
                    hold off; grid on;
                    axis( [ 0, 255, 0, 10*ceil( max( y/10 ) ) ] );
                    pause;
                    close all
                end
                end
                clear x y ym xm km Dx n K xk pk
            else       
                disp 'Empty event!'
                Mu(j,i)=0;
                Sigma(j,i)=0;
                Time(j,i)=time((j-1)*time_buff+1);
            end
        end;
        
        fclose('all'); 
        clear time data Evt
        
end;

BgStruct.Ant=ant;
BgStruct.Loop=nbloop;
BgStruct.Time=Time;
BgStruct.Mu=Mu;
BgStruct.Sigma=Sigma;
if FFTIN
  BgStruct.Spectrum=Spectrum; 
  BgStruct.freq=freq;
end

filename=sprintf('Bgdst_%d.mat',nrun);
filename=[BACKDST_PATH filename];
disp(sprintf('BackgroundDSTBuilder: analysis completed. \nNow saving results to %s',filename))
save(filename,'BgStruct');

if CC==1
    cmd = sprintf('rm -rf /sps/trend/R%d',nrun);
    unix(cmd)
end

            
            
        
