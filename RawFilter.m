function [ Struct] = RawFilter( Struct, id, event )
% RAWFILTER Raw signal filter.
%   1st the stationary noise characteristics are extracted from the data sample 
%   assuming little to no signal. Then a threshold analysis is used in order to  
%   flag significantly over threshold events, in diagreement with a Gaussian 
%   noise hypothesis. These candidate samples are further agglomerated in order
%   to produce time windows for signal events within the sample.


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1-a) Unpack settings
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RunSetup = Struct.Setup;
nrun     = RunSetup.Run;

SharedGlobals

if exist( 'id' )
  DISPLAY = 1; 
end
%DISPLAY

NbCoinc      = RunSetup.TotalCoinc;
Detectors    = [ RunSetup.Det.Name ];
DetectorType = [ RunSetup.Det.isScint ];
CoincStruct  = Struct.Coinc;
if ~isfield( Struct, 'rawfilter' )
  Struct.rawfilter.settings.autocreate = 1;
end
settings = Struct.rawfilter.settings;

MultRef=[Struct.Coinc.Mult];

CaracEvt=cell(NbCoinc,length(Detectors),8);
SigmaSave=zeros(NbCoinc,length(Detectors));
MuSave=zeros(NbCoinc,length(Detectors));
Sat=zeros(NbCoinc,length(Detectors));
minBrut = -ones(NbCoinc,length(Detectors));
maxBrut = -ones(NbCoinc,length(Detectors));
NoBox=zeros(NbCoinc,length(Detectors));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 1-b) Set default settings, whenever
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if ~isfield( settings, 'threshold' )
%   settings.threshold       = 4.0;
% end
% if ~isfield( settings, 'granularity' )
%   settings.granularity     = 25e-9; % s
% end
% if ~isfield( settings, 'max_multiplicty' )
%   settings.max_multiplicty = 6;
% end
% if ~isfield( settings, 'max_total_ToT' )
%   settings.max_total_ToT   = 500e-9; % s
% end
% if ~isfield( settings, 'max_block_ToT' )
%   settings.max_block_ToT   = 300e-9; % s
% end
% if ~isfield( settings, 'bounce_ratio' )
%   settings.bounce_ratio    = 0.8;
% end
% if ~isfield( settings, 'inhibit_window ' )
%   settings.inhibit_window  = 200e-9; % s
% end
% if ~isfield( settings, 'noise_mode' )
%   settings.noise_mode      = 1;
% end

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


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 2-a) Loop on detectors
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i=1:length( Detectors )

  if exist( 'id' ) & ( Detectors( i ) ~= id )
    continue
  end
    
%  if DetectorType( i ) == 0 % antennas
    ib = ibuff;
%   else % scintillators
%     continue % No signal analysis for scintillators
%   end
  
  % Time vector
  %===
  tmu = [ 1:ib ]/FSAMPLING*1e6; % micro-seconds
    
  % Open data file
  %===
  fd = OpenFileData( nrun, Detectors( i ) );
     
  % Select tagged events ( part of coincidences )
  %===
  ind_tag = find( CoincStruct.Det.Tag( :, i ) == 1 );
  evt     = CoincStruct.Det.Evt( ind_tag, i );
  if length( evt ) == 0  % No tagged events for this detector
    disp( sprintf( 'No events for detector %d.', Detectors( i ) ) );
    continue
  end
  devt = [ evt( 1 ) diff( evt )' ]';
  disp( sprintf( 'Starting treatment for detector %d (%d events)...', ...
  Detectors( i ), length( ind_tag ) ) );
    

  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % 2-b) Loop on raw events
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ngood = 0;
  tic;
  for j = 1:length( ind_tag )
    
    % Verbose analysis progress
    %===
    decim = floor( length( ind_tag )/10 );
    if decim>0
		if floor( j/decim ) == j/decim
      	disp(sprintf('Done at %2.0f percent.',j/length(ind_tag)*100))
    	end
	end
    
    % Single event mode?
    %===
    if exist( 'event' ) & event ~= evt(j)
      continue
    end
 
    % Get the selected data
    %===
    fseek( fd, ib*( devt( j )-1 ), 'cof' );
    if exist('event')
        fseek( fd, ib*(evt(j)-1),'bof');
    end;
    DataEvt = double( fread( fd, ib, 'uint8' ) );
    %figure,plot(DataEvt)    
    
    if length( DataEvt ) ~= ib
      disp( sprintf( [ ...
        'Error!\n' ...
        '\tData size for event %d on detector %d is %d samples\n' ... 
        '\tShould have been %d\n' ...
        '\tSkipping it.' ], ...
        evt( j ), Detectors( i ), length( DataEvt ), ib ) );
      continue
    end
    maxBrut(ind_tag(j),i) = max(DataEvt);
    minBrut(ind_tag(j),i) = min(DataEvt);
    if maxBrut(ind_tag(j),i)==255 | minBrut(ind_tag(j),i)==0
        Sat(ind_tag(j),i)=1;
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
    
    if DetectorType(i)==1 % if scintillator
        if sig>1
            sig=0.5;
        end;
    end;
    %mu
    %sig
    MuSave(ind_tag(j),i)=mu;
    SigmaSave(ind_tag(j),i)=sig;
        
    % Window the signal block(s)
    %===
    if DetectorType(i)==0
        
        Win = ( abs( DataEvt - mu ) >= settings.threshold*sig )';

%         modify_flag=zeros(1,length(Win));
%         for k=2:length(Win)
%             if Win(k)==1 & Win(k-1)==0
%                 Win(k-1)=1;
%             elseif Win(k)==0 & Win(k-1)==1 & modify_flag(k-1)==0
%                 Win(k)=1;
%                 modify_flag(k)=1;
%             end;
%         end;
        length(find(Win));
        if sum(Win)>1  % More than 1 point above threshold
          Win = Agglomerate( Win, granularity_sampling);  % join isolated point
        end
        length(find(Win));
    else   % Scintillator       
        Win = ( abs( DataEvt - mu ) >= 2*settings.threshold*sig )';
        Win = Agglomerate( Win, granularity_sampling); 
    end;
    %length(find(Win>0))
    %length(find(Win>0))
    % Get start and end of signal blocks
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
        %DataEvt(iaft:iaft+1*granularity_box)-mu
        %(settings.threshold-2)*sig
        sel = find(abs(DataEvt(iaft:min(length(DataEvt),iaft+1*granularity_box))-mu) >= (settings.threshold-2)*sig);        
        selind = iaft+sel-1;
        selind = selind(selind<=ib);
        Win(selind) = 1;
        %iaft+sel-1
        %pause
      end
      %block_start = max(1,S( bg )-1);
      %block_end   = min(ibuff,S( ed )+1) ;
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
%     if DetectorType(i)==0
%         K = ( block_len >= granularity );
%         n_blocks = sum( K );
%         block_start = block_start( K );
%         block_end   = block_end( K );
%         block_len   = block_len( K );
%     end;
    
    %n_blocks
    %pause;
    %close all
    
    % Global cuts
    %===
    time_over_threshold = sum( block_len )/FSAMPLING;
    
    %%%%%%% On enleve les cuts pour une analyse a posteriori
    if n_blocks == 0 %( n_blocks > settings.max_multiplicty | time_over_threshold > settings.max_total_ToT )
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
    
%       K = ( block_amp >= settings.bounce_ratio*block_amp( 1 ) ) & ...
%           ( block_len <= settings.max_block_ToT ) & ...
%           ( block_dt >= settings.inhibit_window );
%     
%       n_blocks = sum( K );
%       block_start = block_start( K );
%       block_end   = block_end( K );
    end
    
    % Update statistics
    %===
    if n_blocks >= 1
      ngood = ngood + 1;
      is_selected = 1;
    else
      is_selected = 0;
      MultRef(ind_tag(j))=MultRef(ind_tag(j))-1;
      NoBox(ind_tag(j),i)=1;
    end

    CellEvt={is_selected n_blocks time_over_threshold block_start block_end block_len block_amp block_dt};
    
    CaracEvt(ind_tag(j),i,:)=CellEvt;
    
    clear CellEvt
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 2-c) Plot the signal
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %DISPLAY
    if DISPLAY
      disp 'DISPLAY'
      % Title
      %===
      figure( 12 );
      if DetectorType( i ) == 0
        set( 12, ...
          'Name', sprintf( 'Antenna %d - Event %d ', Detectors( i ), evt( j ) ), ...
          'NumberTitle', 'off' );
      elseif DetectorType( i ) == 1
        set( 12, ...
          'Name', sprintf( 'Scintillator %d - Event %d ', Detectors( i ), evt( j ) ), ...
          'NumberTitle', 'off' );
      end
      
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
      pause
    end
  end
  
  % Finalise
  %===
%   dt = toc;
%   fprintf( 'Treatment completed for detector %d.\n', Detectors( i ) );
%   fprintf( '------------------------------------\n' );
%   fprintf( '| Number of triggers  | %10d |\n',   length( ind_tag ) );
%   fprintf( '| Time/event [ mus ]  | %10.0f |\n', 1e6*dt/length( ind_tag ) );
%   fprintf( '| Efficiency [ %% ]    | %10.1f |\n', 100*ngood/length( ind_tag ) );
%   fprintf( '------------------------------------\n' );
  fclose( fd );
end

% Rejected coincidence
Reject=zeros(length(MultRef),1);
ind=find(MultRef<4);
Reject(ind)=1;


Struct.Coinc.Reject.RawFilter=Reject;
Struct.Coinc.Reject.NewMult=MultRef;
Struct.Coinc.Reject.CaracEvt=CaracEvt;
Struct.Coinc.Reject.NoBox=NoBox;

Struct.Coinc.Det.Sigma=SigmaSave;
Struct.Coinc.Det.Mu=MuSave;
Struct.Coinc.Det.MinRaw = minBrut;
Struct.Coinc.Det.MaxRaw = maxBrut;
Struct.Coinc.Det.Sat=Sat;
  
