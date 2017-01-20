function [] =LoopDelay()
% Jianli 17/09/2012

SharedGlobals;

num=[3659 3600 3598 3588 3584 3043 3041 3005 2561];
num = fliplr(num);
AllDetectors=[101:138 140 148:158];
numAntdelay = zeros(length(AllDetectors),length(num));
numAntdelayerr = zeros(length(AllDetectors),length(num));
for i =1 :length(num)
    disp(sprintf('Analysing R%d',num(i))) 
    Antdelay=myPlaneTrackDelayComputation(num(i));
    numAntdelay(:,i)=Antdelay(:,1);
    numAntdelayerr(:,i)=Antdelay(:,2);
    clf;
end

for j=1 : length(AllDetectors) 
    sel = find(abs(numAntdelay(j,:))>0 & abs(numAntdelay(j,:))<50 & abs(numAntdelayerr(j,:))<50);
    if length(sel)>0
        for k = 1:length(sel)
            albl{ k } = num2str( num(k) );
        end
        %
        d = numAntdelay(j,sel);
        e =  numAntdelayerr(j,sel);
        figure(AllDetectors(j))
        errorbar( 1:length(sel), d,e, 'ks','MarkerFaceColor','k','MarkerSize',8 );
        w = 1./e.^2;
        m = mean(d);
        wm = wmean(numAntdelay(j,sel),w);
        wv = var(d,w);
        grid on
        line([1 length(sel)],[m m],'LineWidth',2,'Color','k')
        line([1 length(sel)],[wm wm],'LineWidth',2,'Color','r')
        xlabel('Run ID', labelOpts{:})
        set( gca, 'XTick', 1:length(sel), 'XTickLabel', albl );
        ylabel('Average delay shift [samples]', labelOpts{:})
        disp(sprintf('Antenna %d: mean delay shift = %3.1f pm %3.1f samples',AllDetectors(j),wm,sqrt(wv)))
        text(2,1,sprintf('Antenna %d',AllDetectors(j)))
        saveas(AllDetectors(j),sprintf('Delay_A%d.jpg',AllDetectors(j)),'jpg')
        pause
        close(AllDetectors(j))
    end
end
            
