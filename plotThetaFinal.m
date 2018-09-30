function plotThetaFinal()

 SharedGlobals
 datain = 0;
 if datain
   files = {"thetaplotpi.txt"};
 else
%   files = {"thetaplotpi.txt","thetaplotp.txt","thetaplotpi10.txt"};
   files = {"thetaplotpi.txt"};
 end
 z = rgb('Olive');
 linsty = {'r-','r:','b:'};
 
 for i = 1:length(files)
     a = load(files{i});
     phi = a(:,1);
     N = a(:,2);
     dN = sqrt(a(:,3));

     figure(2)
     errorbar(phi,N,dN, linsty{i},'LineWidth',2 );
     if i == 2
       errorbar(phi,N,dN,':','Color',z,'LineWidth',2 );
     end
     xlabel('Zenith angle (deg)', labelOpts{:})
     ylabel('dN/d\theta (10deg^{-1})', labelOpts{:})
     grid on
     hold on
     xlim([0 90])
 end
 
 if 1
     b = load("candidateswithpsd.txt");
     runs = b(:,2);
     theta = b(:,3);
     phi = b(:,4);
     %plot(phi,40,'k');
     sel = find(runs>3000);
     plotTheta(theta(sel),10,'k');
     %plot(phi,20,'m');
 end
 grid on
 
function plotTheta(theta,step,col)
    SharedGlobals;
    for i=1:90./step
        minedge = (i-1)*step;
        maxedge = i*step;
        sel = find(theta>=minedge & theta<maxedge);
        %pht(i) = mean(phtot(sel));
        tht(i) = (maxedge+minedge)/2;
        N(i) = length(sel);
    end
    figure(2)
    N
    errorbar(tht,N,sqrt(N), 'k-','Color',col,'MarkerFaceColor',col,'LineWidth',2 );
    xlabel('Zenith angle (deg)', labelOpts{:})
    ylabel('Entries (10 deg^{-1})', labelOpts{:})
    grid on
    hold on
    xlim([0 90])