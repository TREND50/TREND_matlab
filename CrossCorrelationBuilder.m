function [ Struct ] = CrossCorrelationBuilder( Struct, id, event )
% Determne trigger time through intercorrelation treatment
% OMH 03/03/11


RunSetup = Struct.Setup;
nrun     = RunSetup.Run;

SharedGlobals;
if exist( 'id' )
  DISPLAY = 1; 
end

NbCoinc      = RunSetup.TotalCoinc;
Detectors    = [ RunSetup.Det.Name ];
DetectorType = [ RunSetup.Det.isScint ];
CoincStruct  = Struct.Coinc;
Evt=[CoincStruct.Det.Evt];
TrigTime = [CoincStruct.Det.TrigTime];
Tag=[Struct.Coinc.Det.Tag];

TrigCor = zeros(NbCoinc,length(Detectors));
CoefCor = zeros(NbCoinc,length(Detectors));

DeltaMin = 100;
DeltaMax = 100;

trig = ibuff/2;  % To be modified if we consider the signal is present at a position different from trigger.
K = trig-DeltaMin:trig+DeltaMax;  % time window centered on trigger time
ts = [ 1: ibuff  ]/FSAMPLING;

for i=1:NbCoinc
  if i/100==floor(i/100)
     disp(sprintf('Coinc %d/%d treated.',i,NbCoinc))
  end
  det_in = find(Tag(i,:)==1 & DetectorType==0);  % Antennas in
  Det = Detectors(det_in);
  EvtNo = Evt(i,det_in);
  TrigRaw = TrigTime(i,det_in);
  L = length(det_in);
  v=zeros(L,ibuff);
  vpeak = zeros(L,1);
  kcorel = zeros(L,1);
  P = zeros(L,1);
  dt=zeros(L*(L-1)/2,1);
  A=zeros(L*(L-1)/2,L);
  l = 1;
  %
  if DISPLAY
      disp(sprintf('Coinc %d: %d antennas',i,L));
  end
  if L<4  % Skip coincs with too few antennas
     continue  
  end
  for j = 1:L-1  % Loop on 1st antenna
    %
    
    fd = OpenFileData( nrun, Det(j));
    if fd==-1
        disp(sprintf('Could not find data file for antenna %d.',Det(j)))
        continue
    end
    fseek(fd,ibuff*(EvtNo(j)-1),'bof');  % A little bit faster
    DataEvt=double(fread(fd,ibuff,'uint8'));
    DataEvt = PassBand( DataEvt, ts, 80.e6, 100.e6 );
    if length(DataEvt)~=ibuff
       disp(sprintf('Error! Data size for event %d on detector %d is %d samples (should be %d)... skipping it.',EvtNo(j),Det(j),length(DataEvt),ibuff));
       break
    end
    fclose(fd);
    if vpeak( j ) == 0
        vpeak( j ) = 0.5*( max(DataEvt) - min(DataEvt) );
        v(j,:) = DataEvt;
    end
      if DISPLAY
          disp(sprintf('\nIntercorrelation with antenna %d\n',Det(j))); 
      end;
      Fint = 10;  % Frequence de sur?chantillonage
      for k = j:L  % Second antenna  
        NK = DeltaMin+DeltaMax+1;
        if vpeak( k ) == 0
            fd = OpenFileData( nrun, Det(k) );
            if fd<0
                disp(sprintf('Could not find data file for antenna %d',Det(k)))
                continue
            end
            fseek(fd,ibuff*(EvtNo(k)-1),'bof');  % A little bit faster
            DataEvt=double(fread(fd,ibuff,'uint8'));
            if length(DataEvt)~=ibuff
              disp(sprintf('Error! Data size for event %d on detector %d is %d samples (should be %d)... skipping it.',EvtNo(k),Det(k),length(DataEvt),ibuff));
              k = j+1;
              break
            end
            fclose(fd);
            v(k,:) = DataEvt;
            vpeak( k ) = 0.5*( max(DataEvt) - min(DataEvt) );
        end
        tmp = ( v(k,K) - mean( v(k,K) ) )/vpeak( k );
        % Interpolation: values of signal vc at points 1:1/Fint:NK
        if ~isfinite(max(tmp))
            disp(sprintf('Warning, no data for antenna %d',Det(k)));
            continue;
        end
        tmp = interp1( 1:NK, tmp, 1:1/Fint:NK, 'spline' );
        dK_filt{ k } = abs( tmp );
        P( k )  = sum( dK_filt{ k }.^2 );

        if k>j
            NK = length( dK_filt{ j } );
            tdiff = ( [ 1:2*NK-1 ]' - NK )/Fint;
            % Intercorrelation between signal l and k
            xc12  = xcorr( dK_filt{ k }, dK_filt{ j } )/sqrt( P( j )*P( k ) );
            % Best correction
            [ cmax, kmax ] = max( abs( xc12 ) );
            %Jx = find( abs( xc12 ) >= 0.99*cmax );
            %kmax = Jx( 1 );
            tmax = tdiff( kmax );
            cmax = xc12(  kmax );

            %% Plots
            if DISPLAY
                figure( 2)
                set(2,'Name',sprintf('Antennas %d & %d ',Det(j),Det(k))); 

                % Raw signals
                subplot( 3, 1, 1 );
                plot( v(j,:), 'blue' );
                hold on;
                plot( v(k,:), 'green' );
                hold off
                grid on


                figure( 2 ); 
                subplot( 3, 1, 2 );
                plot( [ 1:NK ]/Fint, dK_filt{ j }, 'k--','MarkerSize', 2 );
                hold on
                % Shift for plot
                tmaxp=round(tmax*10);
                dKf=dK_filt{ k };
                avi = mean(dK_filt{ k });
                dk_plot =  [ ones( 1, -min( tmaxp, 0 ) )*avi, dKf( max( tmaxp+1, 1 ):min( NK+tmaxp, NK ) ), ones( 1, max( tmaxp, 0 ) )*avi ];  
                plot( [ 1:NK ]/Fint, dk_plot, 'r--','MarkerSize', 2 );
                %plot( [ 1:NK ]/Fint, dK_filt{ k }, 'g--','MarkerSize', 2 );
                hold off; 
                grid on;

                figure( 2 ); 
                subplot( 3, 1, 3 );
                plot( tdiff/FSAMPLING*C0, xc12, 'k' ); 
                hold on;
                plot( tmax/FSAMPLING*C0, cmax, 'ro', 'MarkerFace', 'r', 'MarkerSize', 6 );
                hold off; 
                grid on;

                disp(sprintf( 'coinc %3d: Antenna %4d event %4d | Antenna %4d event %4d\n', i, Det( j ), EvtNo( j ), Det( k ), EvtNo( k ) ));
                disp(sprintf( 'Correction to apply to original trigger times: deltaT   = %6.1f m = %6.1f samples\n', tmax/FSAMPLING*C0, tmax ));
                disp(sprintf( 'corr = %12.1f %%\n', 100*cmax ));
                
                pause

            end
            % Store corrections in vector
            dt(l)=tmax;
            A(l,j)=-1;
            A(l,k)=+1;
            kcorel(j)=kcorel(j)+cmax;
            kcorel(k)=kcorel(k)+cmax;
        end  
        l = l+1;
      end  % end for this pair
  end  % end for this antenna

  if sum(dt)~=0  % This coinc has been analysed
      B = pinv( A )*dt;  
      trigCor = TrigRaw'+B; % Trigger times corrected from cable delays and intercorellation treatment
      trigCor = trigCor-min(trigCor); % Reference is now first antnna truelly triggering: trrigCor is a vector with delays in respect to the 1rst antenna
      TrigCor(i,det_in) = trigCor;
      kcorel=kcorel/(L-1);
      CoefCor(i,det_in) = kcorel;
  end
end  % End for this coinc
Struct.Coinc.Det.TrigCor = TrigCor;
Struct.Coinc.Det.CoefCor = CoefCor;
