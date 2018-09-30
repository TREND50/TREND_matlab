function plotPhi()

 a = load("candidates.txt");
 phi = a(:,4);
 plot(phi,40,'b');
 %plot(phi,20,'b');
 runs = a(:,2);
 sel = find(runs>3000);
 plot(phi(sel),20,'r');
 disp "Nb events with",length(sel)
 
 b = load("candidateswithpsd.txt");
 runs = b(:,2);
 phi = b(:,4);
 plot(phi,40,'k');
 sel = find(runs>3000);
 plot(phi(sel),40,'g');
 %plot(phi,20,'m');
 

function plot(phtot,step,col)
    SharedGlobals;
    for i=1:360./step
        minedge = (i-1)*step;
        maxedge = i*step;
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
    N = N*20/step;
    figure(1)
    errorbar(pht,N,sqrt(N), 'r-','Color',col,'MarkerFaceColor',col,'LineWidth',2 );
    xlabel('Azimuth angle [deg]', labelOpts{:})
    ylabel('dN/d\phi [20deg^{-1}]', labelOpts{:})
    grid on
    hold on
    xlim([0 360])