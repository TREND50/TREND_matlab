function [CaracEvt] = PostRawFilter( Struct, id )
% RAWFILTER Raw signal filter for pulses from Candidates
% Taken from RawFilter.m
% Treatment for 
% OMH 04/12/2012

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1-a) Unpack settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RunSetup = Struct.Setup;
nrun     = RunSetup.Run;

SharedGlobals

if exist( 'id' )
  DISPLAY = 0; 
end

NbCoinc      = RunSetup.TotalCoinc;
Detectors    = [ RunSetup.Det.Name ];
DetectorType = [ RunSetup.Det.isScint ];
CoincStruct  = Struct.Coinc;
if ~isfield( Struct, 'rawfilter' )
  Struct.rawfilter.settings.autocreate = 1;
end
settings = Struct.rawfilter.settings;

MultRef=[Struct.Coinc.Mult];

CaracEvt=cell(length(Detectors),8);
SigmaSave=zeros(1,length(Detectors));
MuSave=zeros(1,length(Detectors));
Sat=zeros(1,length(Detectors));
minBrut = -ones(1,length(Detectors));
maxBrut = -ones(1,length(Detectors));
NoBox=zeros(1,length(Detectors));

settings.threshold = threshold;
settings.granularity_box = granularity_box; % s
settings.granularity_sampling = granularity_sampling; % s
settings.max_out = max_out;
settings.max_total_ToT   = max_out_ToT; % s
settings.max_block_ToT   = max_block_ToT; % s
%settings.bounce_ratio    = bounce_ratio;
%settings.inhibit_window  = inhibit_window; % s
settings.noise_mode      = noise_mode;

Struct.rawfilter.settings = settings; % update structure
granularity_box = round( settings.granularity_box*FSAMPLING );
granularity_sampling = round( settings.granularity_sampling*FSAMPLING );

Evt = CoincStruct.Det.Evt;
signals = CoincStruct.ShowerRecons.ShowerSignals;
CoincId = CoincStruct.IdCoinc;
ind = find(CoincId==id);
trig = CoincStruct.Det.Tag(ind,:);
nTrig = sum(trig);
indant = find(trig==1);

ib = ibuff;
tmu = [ 1:ib ]/FSAMPLING*1e6; % micro-seconds
ngood = 0;
    
%% Loop on detectors
for i = 1:nTrig    

    DataEvt = signals{ind,indant(i)}/SCALE;  % Amplitude back in LSB
    
    %figure,plot(DataEvt)    
    
    if length( DataEvt ) ~= ib
      disp( sprintf( [ ...
        'Error!\n' ...
        '\tData size for event %d on detector %d is %d samples\n' ... 
        '\tShould have been %d\n' ...
        '\tSkipping it.' ], ...
        Evt(ind,indant(i)), Detectors( indant(i) ), length( DataEvt ), ib ) );
      continue
    end
    maxBrut(indant(i)) = max(DataEvt);
    minBrut(indant(i)) = min(DataEvt);
    if maxBrut(indant(i))==255 | minBrut(indant(i))==0
        Sat(indant(i))=1;
    end;
        
    % Extract the stationary noise parameters
    %===
    KLOW = [ 1:(round( length( DataEvt )/2 )-100) ]; % 1st half subset of data. Before trigger
    indmax=find(max(DataEvt(KLOW))==DataEvt(KLOW));
    if DataEvt(KLOW(indmax))>155 & indmax>150
        KLOW=[1:indmax-10];
    end;
    if ( settings.noise_mode == 1 ) | ( DISPLAY == 1 )
      [ y, x ] = hist( DataEvt( KLOW ), [ 0:255 ] );
      [ ym, km ] = max( y );
      xm  = x( km );
      Dx  = unique( abs( x( y > 0.5*ym ) - xm ) ); % Select samples within Width at Half Full Height (WHFH).
      n = round( Dx(end)*2.5577 ); % 3 sigmas width, corrected from WHFH and rounded.  
      K   = km + [ -n:1:n ];
      K   = K( K <= 256 & K >= 1 ); 
      xk  = x( K );
      pk  = y( K )/sum( y( K ) ); 
      mu  = sum( xk.*pk )  ;               % Noise mean value
      if size(Dx,2)<4  % Amplitude histogram too sharp
        sig=std(DataEvt(KLOW));
      else
        sig = sqrt( sum( ( xk - mu ).^2.*pk ) ); % Noise sigma
      end
    end
    
    if ( settings.noise_mode ~= 1 )
      mu  = mean( DataEvt( KLOW ) );
      sig = std( DataEvt( KLOW ) );
    end
    MuSave(indant(i))=mu;
    SigmaSave(indant(i))=sig;
        
    % Window the signal block(s)
    %===
    Win = ( abs( DataEvt - mu ) >= settings.threshold*sig )';
    % Get start    Win = ( abs( DataEvt - mu ) >= settings.threshold*sig )';
    if sum(Win)>1  % More than 1 point above threshold
      Win = Agglomerate( Win, granularity_sampling);  % join isolated point
    end

    %===
    S  = find( Win );
    if isempty( S )
      n_blocks = 0;
      block_start = [];
      block_end   = []; 
    else
      dS = [ 2, diff( S ) ];
      bg = find( dS > 1 );
      ed = [ bg(2:end)-1, length( S ) ];
  
      n_blocks    = length( bg );
      block_start = S(bg);
      block_end = S(ed);
      for k=1:n_blocks
        %k
        ibef = block_start(k)-1;
        %DataEvt(ibef-1*granularity_box:ibef)-mu
        %DataEvt(sel+block_start(k)-granularity_sampling)-mu
        %(settings.threshold-2)*sig
        sel = find(abs(DataEvt(max(1,ibef-granularity_box):ibef)-mu) >= (settings.threshold-2)*sig);
        selind = ibef-1*granularity_box+sel-1;
        selind = selind(selind>0);
        Win(selind) = 1;
        %pause
        iaft = block_end(k)+1;
        sel = find(abs(DataEvt(iaft:min(length(DataEvt),iaft+1*granularity_box))-mu) >= (settings.threshold-2)*sig);        
        selind = iaft+sel-1;
        selind = selind(selind<=ib);
        Win(selind) = 1;
        %iaft+sel-1
        %pause
      end
    end
    
    %Win = Agglomerate( Win, granularity_sampling);  % join isolated point
    Win = Agglomerate( Win, granularity_box);   % join box
    S  = find( Win );
    if isempty( S )
      n_blocks = 0;
      block_start = [];
      block_end   = []; 
    else
      dS = [ 2, diff( S ) ];
      bg = find( dS > 1 );
      ed = [ bg(2:end)-1, length( S ) ];
      n_blocks    = length( bg );
      block_start = S(bg);
      block_end = S(ed);
    end
        
    % Clean blocks shorter than granularity
    %===
    block_len = block_end - block_start + 1;
    
    % Global cuts
    %===
    time_over_threshold = sum( block_len )/FSAMPLING;
    
    %%%%%%% On enleve les cuts pour une analyse a posteriori
    if n_blocks == 0 
      n_blocks    = 0;
      block_start = [];
      block_end   = [];
      block_len   = [];
      block_amp   = [];
      block_dt    = [];
    elseif ( n_blocks > 0 )   
      % Clean reflections and blocks obviously too long
      %===
      block_len = block_len/FSAMPLING;
      block_amp = zeros( size( block_len ) );
      for k = 1:n_blocks
        K = [ block_start( k ):1:block_end( k ) ];
        block_amp( k ) = max( abs( DataEvt( K ) - mu ) );
      end
      
      if n_blocks >= 2
        block_dt = [ 1.0, ( block_start( 2:end ) - block_end( 1:end-1 ) )/FSAMPLING ];
      else
        block_dt = [ 1.0 ];
      end
    end
    
    % Update statistics
    %===
    if n_blocks >= 1
      ngood = ngood + 1;
      is_selected = 1;
    else
      is_selected = 0;
      MultRef(indant(i))=MultRef(indant(i))-1;
      NoBox(indant(i))=1;
    end

    CellEvt={is_selected n_blocks time_over_threshold block_start block_end block_len block_amp block_dt};
    
    CaracEvt(indant(i),:)=CellEvt;
    
    clear CellEvt
    
    %% DISPLAY
    if DISPLAY
      disp 'DISPLAY'
      figure( 12 );
      set( 12,'Name', sprintf( 'Event %d - Antenna %d ', Evt(ind,indant(i)), Detectors(indant(i))),'NumberTitle', 'off' );
   
      % ADCs distribution
      %===
      subplot( 2, 1, 1 );
      h = bar( x, y, 1.0 ); hold on;
      set( h, 'FaceColor', 'none', 'EdgeColor', 'k' );
      plot( x, ym*exp( -0.5*( x - mu ).^2/sig^2 ), 'r', 'LineWidth', 2 );
      hold off; grid on;
      axis( [ 0, 255, 0, 10*ceil( max( y/10 ) ) ] );
      xlabel( 'ADCs [ LSB ]', labelOpts{:} );
      ylabel( '# entries', labelOpts{:} );
      set( gca, 'FontSize', 14 );
      strtitle = { 'REJECT', 'ACCEPT' };
      title( strtitle{ is_selected+1 }, 'FontSize', 20 );
      
      % Signal as time
      %===
      subplot( 2, 1, 2 );
      tmu = 1:length(DataEvt);
      plot( tmu, DataEvt, 'k-' ); hold on;
      for k = 1:n_blocks
        t1 = tmu( block_start( k ) );
        t2 = tmu( block_end( k ) );
        %t1 =  block_start( k );
        %t2 =  block_end( k );
        lineOpts = { 'Color', 'r', 'LineWidth', 1 };
        line( t1*[ 1, 1 ], [ 0, 255 ], lineOpts{:} );
        line( t2*[ 1, 1 ], [ 0, 255 ], lineOpts{:} );
        line( [ t1, t2 ], [ 0, 0 ], lineOpts{:} );
        line( [ t1, t2 ], [ 255, 255 ], lineOpts{:} );
      end
      hold off; grid on;
      axis( [ min( tmu ), max( tmu ), 0, 255 ] );
      line( [ min(tmu),max(tmu) ], [ mu + settings.threshold*sig mu + settings.threshold*sig ], 'Color','b' );
      line( [ min(tmu),max(tmu) ], [ mu - settings.threshold*sig mu - settings.threshold*sig ], 'Color','b' );
      line( [ min(tmu),max(tmu) ], [ mu + (settings.threshold-2)*sig mu + (settings.threshold-2)*sig ], 'Color','g' );
      line( [ min(tmu),max(tmu) ], [ mu - (settings.threshold-2)*sig mu - (settings.threshold-2)*sig ], 'Color','g' );
      line( [ min(tmu),max(tmu) ], [ mu + (settings.threshold+1)*sig mu + (settings.threshold+1)*sig ], 'Color','y' );
      line( [ min(tmu),max(tmu) ], [ mu - (settings.threshold+1)*sig mu - (settings.threshold+1)*sig ], 'Color','y' );
      
      xlabel( 'Time [ \mus ]', labelOpts{:} );
      ylabel( 'ADCs [ LSB ]', labelOpts{:} );
      set( gca, 'FontSize', 14 );
          
      sig
      n_blocks 
      ToT = time_over_threshold*1e9 
      Dur = block_len*1e9 
      block_amp
      CaracEvt
      pause
    end
    
end

  
