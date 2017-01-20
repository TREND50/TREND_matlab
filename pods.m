function []=pods(nrun,antennas,R)
% Function to plot pods position
% + antennas layout for run nrun if nrun exists.
% OMH 17/09/09


if nrun == 2000
    shift = 0;  % Shift all antennas North?
end

%% Plot options
%%===
labelOpts = { 'FontSize', 18 };
axisOpts  = { 'fontSize', 14, 'YScale', 'linear' };
mtyp= ['dblu';'pbla';'hgre'];
ltyp= ['blu';'bla';'gre'];
scrsz = get( 0, 'ScreenSize' );

[ nevent, init_ant ] = get_nevent( nrun,0 );
%init_ant = 135;

%% Load pods position
pods = load('pods.txt');
pods_n = pods(:,1);
pods_we = pods(:,2);
pods_sn = pods(:,3);
pods_alt = pods(:,4);
pods_xp = pods(:,5);

%% Check altitudes
% delta = pods_alt-pods_xp;
% figure(12)
% subplot(1,2,1)
% hist(delta)
% xlabel('DeltaAlt [m]',labelOpts{:})
% subplot(1,2,2)
% plot(pods_n,delta,'*')
% grid on
% line([120.5 120.5],[-10 10])
% line([140.5 140.5],[-10 10])
% line([160.5 160.5],[-10 10])
% xlabel('Pod ID',labelOpts{:})
% ylabel('DeltaAlt [m]',labelOpts{:})
%for i=1:length(pods_we)
%    z=getElevation(pods_we(i),pods_sn(i));
%    disp(sprintf('%d %3.1f',pods_n(i),z));
%end

%% Load antenna positions
if exist('nrun')
    position_ant_all   = open_coordonnees( nrun );
    position_ant = position_ant_all(find(position_ant_all(:,4)>0),:);
else
    position_ant = zeros(1,4);
end
if nrun == 2000
    position_ant(:,3) = position_ant(:,3)+shift;
end

%% Closest neighbours
pos = position_ant(find(position_ant(:,4)>0),2:4);
dis=zeros(size(pos,1),size(pos,1));
for i=1:size(pos,1)
    for j=1:size(pos,1)
        dist(i,j)=norm(pos(i,:)-pos(j,:));
    end
    dist(i,i)=1e3;
    [mdist(i) idist] = min(dist(i,:));
    adist(i) = position_ant(idist,1);
end
% disp(sprintf('%d antennas',length(adist)));
% disp 'Distance to closest neighbour:'
% for i=1:size(position_ant,1)
%     disp(sprintf('A%d-A%d : %3.1f m',position_ant(i,1),adist(i),mdist(i)));
% end
%% Distance to pods
for i=1:size(position_ant,1)
    a(i)=find(pods_n==position_ant(i,1));
end
dist_pod = sqrt((position_ant(:,2)-pods_we(a)).^2+(position_ant(:,3)-pods_sn(a)).^2);
disp 'Distance to connected pod:'
for i=1:size(position_ant,1)
    disp(sprintf('A%d : %3.1f m',position_ant(i,1),dist_pod(i)));
end
figure(12)
hist(dist_pod)
datastats(dist_pod)
disp(sprintf('Total cable length = %3.1f m ', sum(dist_pod) ))
length(find(dist_pod<150))
%% Scintillators
if exist('nrun')&nrun>1727
    position_scint = open_coordonnees_scint( nrun );
    position_scint = position_scint(find(position_scint(:,4)>0),:);
else
    position_scint = [];
end

%% Plot
figure(1)
tit='21CMA';     
tit = '';
ftit=sprintf('21CMA.jpg');
scrsz = get( 0, 'ScreenSize' );
set(1, 'Name', tit,'NumberTitle','off','Position',[1 scrsz(2) scrsz(3)/3 scrsz(4)]);
title(tit,labelOpts{:})
grid on;
hold on;

%% Load relief map
if nrun>1000 & nrun<2000  % CrossPoint
    xmin = -100;
    xmax = 400;
    ymin = -500;
    ymax = 700;
elseif nrun == 1000
    xmin = -100;
    xmax = 400;
    ymin = -500;
    ymax = 800;
elseif nrun == 2000 | nrun>2400 %East arm
    %xmin = -100;
    xmin = -4500;
    %xmax = 3000;
    xmax = 4500;
    %ymin = -500;
    ymin = -2000;
    %ymax = 800;
    ymax = 7000;
    %ymax = 4000;
    set(1, 'Name', tit,'NumberTitle','off','Position',[1 1 scrsz(3) scrsz(4)/2]);
else
    xmin = -100;
    xmax = 3000;
    ymin = -500;
    ymax = 4000;
end
[x,y,z] = getMap('map.bin',[xmin,xmax,ymin,ymax]);
axis([xmin,xmax,ymin,ymax])
    
contourf(x,y,z,40)
colormap bone
colorbar
%axis equal

%% Pods:
px=3.031;
py=2.57;
xpod=[-3*px -3*px 0 3*px 3*px 0 -3*px];
ypod=[3*py -3*py -6*py -3*py 3*py 6*py 3*py];

for i=1:length(pods_n)
    %plot(pods_we(i),pods_sn(i),'dw','MarkerSize',8);
    plot(xpod+pods_we(i),ypod+pods_sn(i),'w');
end
xlabel( 'W-E [ m ]', labelOpts{:} );
ylabel( 'S-N [ m ]', labelOpts{:} );

%if size(position_ant,1)<2
if 1
    for i = 1:length( pods_n )
        if nrun==2000 & pods_we(i) > -10000  % display WE arm only
            %continue
        end
        sign = 1-4*(i/2-floor(i/2));       
        if pods_n(i)<141
            text( pods_we(i), pods_sn(i)-7+sign*50, num2str( pods_n( i ) ), 'FontSize', 8, 'FontWeight', 'bold', 'Color','g' );
        else
            text( pods_we(i)-20+sign*50, pods_sn(i), num2str( pods_n( i ) ), 'FontSize', 8, 'FontWeight', 'bold', 'Color','g' );
        end
    end
end

%% Antennas
if nrun>0
    plot( position_ant( find(dist_pod>0), 2 ), position_ant( find(dist_pod>0), 3 ), '^y', 'MarkerSize', 8, 'MarkerFaceColor','y' );  %Ground view  
    % Warning!!! Different conventions from what was used so far! 
    % 2nd col: WE (was SN)
    % 3rd col: SN  (was EW))
    for i = 1:size( position_ant,1 )        
        text( position_ant( i , 2 )+20, position_ant( i, 3 ), num2str( position_ant( i, 1 ) ), 'FontSize', 12, 'FontWeight', 'bold','Color','y' );
    end
end

%% Scintillators
if size(position_scint,1)>0
    %plot( position_scint( :, 2 ), position_scint( :, 3 ), 'sr', 'MarkerSize', 12, 'MarkerFaceColor','r' );  %Ground view  
    % Warning!!! Different conventions from what was used so far! 
    for i = 1:size( position_scint,1 )
        %text( position_scint( i , 2 )+20, position_scint( i, 3 ), num2str( position_scint( i, 1 ) ), 'FontSize', 12, 'FontWeight', 'bold','Color','k' );
    end
end
if 0
nscint = 0;
xscint = -50:150:550;
yscint = -200:150:100;
for i = 1:length(xscint)
    for j = 1:length(yscint)
        plot(xscint(i),yscint(j),'sr','MarkerFaceColor','r')
        nscint = nscint+1;
    end
end
xscint = -50:100:150;
yscint = 250:150:600;
for i = 1:length(xscint)
    for j = 1:length(yscint)
        plot(xscint(i),yscint(j),'sr','MarkerFaceColor','r')
        nscint = nscint +1;
    end
end
end
if nrun>999
    %crossPointEnvir
end
grid on
saveas(1,ftit);

%% Plot circles
if exist('antennas','var') & length(antennas)>1
    possort = sortrows(position_ant_all,1);
    for i=1:length(antennas)
        ind = antennas(i)-init_ant;
        myCircle(possort( ind, 2 ),possort( ind, 3 ),R(i));
    end
end

%% 3D plot
% figure(2)
% surf(x,y,z)
% axis([-500 800 -1000 800 2600 2800])
% colormap copper
% colorbar
% hold on
% plot3( position_ant( :, 2 ), position_ant( :, 3 ), getElevation(position_ant( :, 2 ), position_ant( :, 3 )),'b*', 'MarkerSize', 12 );  %Ground view  

% for i=1:length(pods_n)
%     plot3(xpod+pods_we(i),ypod+pods_sn(i),getElevation(pods_we(i),pods_sn(i)),'g');
% end
