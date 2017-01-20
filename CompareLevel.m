function [] = CompareLevel(nrun1,nrun2)
% Compare mean PSDs for 2 nruns
% Taken from 'compare_level
% OMH 08/06/2011

SharedGlobals;  
labelOpts = { 'FontSize', 18 };

antenna = [ 101:138, 140, 148:158 ];
N = length(antenna);
level = zeros(1,N);
stdlevel = zeros(1,N);
Fmin=55e6;
Fmax= 95e6;
for k =1:N
    d1 = ReadPSD( nrun1, antenna(k) );
    d2 = ReadPSD( nrun2, antenna(k) );
    if size(d1,1)==0 | size(d2,1)==0
        disp 'skip'
        continue
    end
    f  = d1.f;  %Frequency
    time1 = d1.t;
    time2 = d2.t;
    BW = ( f >= Fmin & f <= Fmax );    
    %psd1 = mean(d1.psd,2);
    %psd2 = mean(d2.psd,2);
    psd1 = d1.psd(:,end); % Last time measurement of 2nd run
    psd2 = d2.psd(:,1);  % First time measurement of 2nd run
    diff = psd2-psd1;
    BW = ( f >= Fmin & f <= Fmax );
    level(k) = mean(diff(BW));
    stdlevel( k ) = std( diff( BW ) )/sqrt( sum( BW ) );
    albl{ k } = num2str( antenna( k ) );
    
    figure(k)
    set(k,'Name',sprintf('Antenna %d - PSD R%d-R%d',antenna(k),nrun1,nrun2),'NumberTitle','off') 
    subplot(2,1,1)
    plot(f,psd1,'k','LineWidth',2)
    hold on; grid on
    plot(f,psd2,'r','LineWidth',2)
    xlabel('Frequency [MHz]',labelOpts{:})
    subplot(2,1,2)
    plot(f,diff,'g','LineWidth',2)
    grid on;
    xlabel('Frequency [MHz]',labelOpts{:})
    %pause
    close(k)
end

ref_level = mean( level( antenna >= 101 & antenna <= 119 ) );
figure( 100 ); 
hold on;
plot( [ 0, N+1 ], ref_level*[ 1, 1 ], 'r', 'LineWidth', 2 );
errorbar( 1:N, level, stdlevel, 'ko' );
xlim([0 N+2])
hold off; 
grid on;
set( gca, 'FontSize', 10, 'FontWeight', 'bold' );
set( gca, 'XTick', 1:N, 'XTickLabel', albl );
ylabel( 'Level [ dB  ref  V / Hz^{1/2} ]' );
xlabel( 'Antenna' );