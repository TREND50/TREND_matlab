function plotPhiFinal()

 SharedGlobals
 %   files = {"phiplotpi.txt","phiplotp.txt","phiplotpi10.txt"};
 files = {"phiplotpi20deg.txt"};
 z = rgb('Olive');
 linsty = {'r-','r:','b:'};
 
 for i = 1:length(files)
     a = load(files{i});
     phi = a(:,1);
     N = a(:,2);
     dN = sqrt(a(:,3));

     figure(1)
     errorbar(phi,N,dN, linsty{i},'LineWidth',2 );
     if i == 2
       errorbar(phi,N,dN,':','Color',z,'LineWidth',2 );
     end
     xlabel('Azimuth angle (deg)', labelOpts{:})
     ylabel('Entries (40 deg^{-1})', labelOpts{:})
     hold on
     xlim([0 360])
 end
 
 if 1
     b = load("candidateswithpsd.txt");
     runs = b(:,2);
     phi = b(:,4);
     sel = find(runs>3000);
     plot(phi(sel),20,'k');
     %
     b = load("candidates.txt");
     runs = b(:,2);
     phi = b(:,4);
     sel = find(runs>0);
     plot(phi(sel),20,'g');
 end
 
 ax1 = gca; % current axes
 ax1_pos = ax1.Position; % position of first axes
 ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','Color','none');
 x2 = -90:1:270;
 y2 = 0*x2;
 line(x2,y2,'Parent',ax2,'Color','w')
 grid on
 
function plot(phtot,step,col)
    SharedGlobals;
    for i=1:360./step
        minedge = (i-1)*step-10;
        maxedge = i*step-10;
        minsel = mod(minedge-90,360);
        maxsel = mod(maxedge-90,360);
        if minsel>maxsel
           sel = find( (phtot>=minsel | phtot<maxsel));
        else
            sel = find(phtot>=minsel & phtot<maxsel);
        end
        %pht(i) = mean(phtot(sel));
        pht(i) = (maxedge+minedge)/2;
        N(i) = length(sel);
    end
    N = N*408./sum(N);
    figure(1)
    errorbar(pht,N,sqrt(N), 'r-','Color',col,'MarkerFaceColor',col,'LineWidth',2 );
    xlabel('Azimuth angle (deg)', labelOpts{:})
    ylabel('dN/d\phi (40deg^{-1})', labelOpts{:})
    hold on
    xlim([0 360])