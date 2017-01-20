function [plotid]=Layout(nrun,dst)
% Display TREND layout for triggered detectirs in this run
% Adaptated from pods.m
% OMH 03/06/11

SharedGlobals;

%% Plot options
%%===
labelOpts = { 'FontSize', 18 };
scrsz = get( 0, 'ScreenSize' );

%% Load dst
if ~exist('dst')
    dstname = [DST_PATH sprintf(dst_filename,nrun,1)];
    dst = load(dstname);
end
DetStruct = dst.Struct.Setup.Det;
Detectors = [DetStruct.Name];
isScint = [DetStruct.isScint];
xpos_ant = [DetStruct.X];
ypos_ant = [DetStruct.Y];
zpos_ant = [DetStruct.Z];
PodNb = [DetStruct.PodNb];
xpos_pod = [DetStruct.PodX];
ypos_pod = [DetStruct.PodY];
zpos_pod = [DetStruct.PodZ];

if 0
%% Closest neighbours
ants = find(isScint==0);
Antennas= Detectors(ants);
pos = [xpos_ant(ants)' ypos_ant(ants)' zpos_ant(ants)'];
dist=zeros(size(pos,1),size(pos,1));
for i=1:length(Antennas)
    for j=1:length(Antennas)
        dist(i,j)=norm(pos(i,:)-pos(j,:));
    end
    dist(i,i)=1e3;
    [mdist(i) idist] = min(dist(i,:));
    adist(i) = Antennas(idist);
end
disp(sprintf('%d antennas',length(Antennas)));
disp 'Distance to closest neighbour:'
for i=1:length(Antennas)
     disp(sprintf('A%d-A%d : %3.1f m',Antennas(i),adist(i),mdist(i)));
end

%% Distance to pods
pos_det = [xpos_ant' ypos_ant' zpos_ant'];
pos_pod = [xpos_pod' ypos_pod' zpos_pod'];
for i = 1:length(Detectors)
    dist_pod(i) = norm(pos_det(i,:)-pos_pod(i,:));
end
disp(sprintf('%d detectors',length(Detectors)));
disp 'Distance to connected pod:'
for i=1:length(Detectors)
    if isScint(i)==1
        disp(sprintf('Scint %d : %3.1f m',Detectors(i),dist_pod(i)));
    else
        disp(sprintf('A%d : %3.1f m',Detectors(i),dist_pod(i)));
    end
end
figure(12)
hist(dist_pod)
datastats(dist_pod)
disp(sprintf('Total cable length = %3.1f m ', sum(dist_pod) ))

end

%% Plot
plotid = nrun;
figure(plotid)
xmin = -100;
%xmin = -4500;
xmax = 3000;
%xmax = 4500;
ymin = -500;
%ymin = -2000;
ymax = 800;
%ymax = 7000;
set(plotid, 'Name', sprintf('TREND layout - Run %d', nrun),'NumberTitle','off','Position',[1 1 scrsz(3)/1.5 scrsz(4)/2]);
axis([xmin,xmax,ymin,ymax])
xlabel( 'W-E [ m ]', labelOpts{:} );
ylabel( 'S-N [ m ]', labelOpts{:} );
hold on    
grid on

%% Load relief map
%[x,y,z] = getMap('map.bin',[xmin,xmax,ymin,ymax]);
%contourf(x,y,z,40)
%colormap bone
%colorbar
%axis equal

%% Pods:
px=3.031;
py=2.57;
xpod=[-3*px -3*px 0 3*px 3*px 0 -3*px];
ypod=[3*py -3*py -6*py -3*py 3*py 6*py 3*py];
for i=1:length(Detectors)
    %plot(xpos_pod(i),ypos_pod(i),'dw','MarkerSize',8);
    plot(xpod+xpos_pod(i),ypod+ypos_pod(i),'w');
    if 0
    sign = 1-4*(i/2-floor(i/2));       
    if PodNb(i)<141
        text( xpos_pod(i), ypos_pod(i)-7+sign*50, num2str( PodNb( i ) ), 'FontSize', 8, 'FontWeight', 'bold', 'Color','g' );
    else
        text( xpos_pod(i)-20+sign*50, ypos_pod(i), num2str( PodNb( i ) ), 'FontSize', 8, 'FontWeight', 'bold', 'Color','g' );
    end
    end
end

%% Detectors
for i = 1:length( Detectors )
    if isScint(i)==1
        plot( xpos_ant(i), ypos_ant(i), 'sr', 'MarkerSize', 8, 'MarkerFaceColor','r' );  %Ground view
        text( xpos_ant(i)+20, ypos_ant(i), num2str( Detectors(i) ), 'FontSize', 12, 'FontWeight', 'bold','Color','r' );
    else
        plot( xpos_ant(i), ypos_ant(i), '^k', 'MarkerSize', 8, 'MarkerFaceColor','y' );  %Ground view
        text( xpos_ant(i)+20, ypos_ant(i), num2str( Detectors(i) ), 'FontSize', 12, 'FontWeight', 'bold','Color','k' );
    end
end
