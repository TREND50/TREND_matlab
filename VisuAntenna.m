function [fov] = VisuAntenna(ant)
% Determine the angular position of horizon as a function of azimuth for antenna ant at position [x,y,z] value of 
% Taken from VisuAntenna by Dimitri
% OMH 03/02/2016

SharedGlobals;
DISPLAY = 2
[pos_ant] = SetupCharacteristics_TREND50(ant,3003);
x = pos_ant(1);
y = pos_ant(2);
z = pos_ant(3);
disp(sprintf('VisuAntenna.m: computing horizon elevation for antenna %d at position (%3.1f, %3.1f, %3.1f)', ant, x, y, z))

fovname = sprintf('AntennaFOV_A%d.mat',ant);
if 0
    %exist(fovname)
   display(sprintf('Loading %s...',fovname)) 
   load(fovname);
   return
end
disp 'Computing FOV...'
%DISPLAY = 0;
[xref,yref,zref] = getMap('map.bin',[-150e3,+150e3,-150e3,+150e3]);
step=100;
dmax=20000;
phi=[1:1:360];

AntennaBlind=zeros(length(ant),length(phi));

%phi=100;

if DISPLAY>0
    dat=figure(1);
    scrsz=get(0,'ScreenSize');
    set(dat,'Position',[scrsz(1) scrsz(2) scrsz(3)*0.99 scrsz(4)*0.85]);
    dat1=subplot(1,2,1);
    plot(x,y,'+')
    xlim([-10000,+10000])
    ylim([-10000,+10000])
    ref=line(0,0,'color','k','Marker','+','erasemode','xor');
    refphi=line(0,0,'color','r','Marker','+','erasemode','xor');
    refmap=line(0,0,'color','r','Marker','o','erasemode','xor');
    dat2=subplot(2,2,2);
    altprofile=line(0,0,'color','b','Marker','+','erasemode','xor');
    altprofilestruct=line(0,0,'color','k','Marker','.','erasemode','xor');
    set(dat2,'XLim',[0 dmax]);
    dat3=subplot(2,4,7);
    set(dat3,'XLim',[0 360],'YLim',[0 20])
    blindangle=line(0,0,'color','b','Marker','+','erasemode','xor');
    blindanglestruct=line(0,0,'color','k','Marker','.','erasemode','xor');
end

% Loop on antennas 
for i=1:length(ant)
    % Loop on phi angle
    for j=1:length(phi)
        if DISPLAY>0
            set(refphi,'xdata',[x(i) x(i)+dmax*sind(-phi(j))],'ydata',[y(i) y(i)+dmax*cosd(-phi(j))])
        end
        % Scan track
        for k=1:(dmax/step)
            
            X(k)=x(i)+k*step*sind(-phi(j));
            Y(k)=y(i)+k*step*cosd(-phi(j));
            %indx=find(abs(X(k)-xref)==min(abs(X(k)-xref)));
            [a indx]=min(abs(X(k)-xref));
            if length(indx>1)
                indx=indx(1);
            end;
            %indy=find(abs(Y(k)-yref)==min(abs(Y(k)-yref)));
            [a indy]=min(abs(Y(k)-yref));
            if length(indy>1)
                indy=indy(1);
            end;
            zstruct(k)=zref(indy,indx);  % Watch out for index order!!! y THEN x!
            dist(k)=sqrt((x(i)-X(k))^2+(y(i)-Y(k))^2);
            Zstruct(k)=zstruct(k)-z;
            ztotal(k)=getElevation(X(k),Y(k));
            Z(k)=ztotal(k)-z;
            %[xref(indx) X(k) yref(indy) Y(k) zref(indx,indy) zref(indy,indx) ztotal(k)]
            %pause
            if Zstruct(k)>0
                %thetablind(k)=atand(Z(k)/dist(k));
                thetablindstruct(k)=atand(Zstruct(k)/dist(k));
            else
                %thetablind(k)=0;
                thetablindstruct(k)=0;
            end;

            if DISPLAY>1
                set(refmap,'xdata',X(k),'ydata',Y(k));
                set(altprofile,'xdata',[dist],'ydata',[Z])
                set(altprofilestruct,'xdata',[dist],'ydata',[Zstruct])
                if k>1
                    set(dat2,'YLim',[1.1*min(Z) 1.1*max(Z)])
                end;
                drawnow
                %pause;
            end

        end;
        
%     indblind=find(Z==max(Z));  
%     Zblind(j)=Z(indblind);
%     %Xblind(j)=X(indblind);
%     %Yblind(j)=Y(indblind);
%     Distblind(j)=dist(indblind);
%     if Z(k)>0
%         Thetablind(j)=atand(Zblind(j)/Distblind(j))
%     else
%         Thetablind(j)=0
%     end;

        %Thetablind(j)=max(thetablind);
        Thetablindstruct(j)=max(thetablindstruct);
        AntennaBlind(i,j)=Thetablindstruct(j);
        if DISPLAY
            set(refmap,'xdata',X(k),'ydata',Y(k));
            %set(altprofile,'xdata',[dist],'ydata',[Z])  %blue
            set(altprofilestruct,'xdata',[dist],'ydata',[Zstruct])  %red
            if k>1
              set(dat2,'YLim',[1.1*min(Zstruct) 1.1*max(Zstruct)])
            end;
            %set(blindangle,'xdata',[phi(1:j)],'ydata',[Thetablind(1:j)])
            set(blindanglestruct,'xdata',[phi(1:j)],'ydata',[Thetablindstruct(1:j)])
            drawnow
            %pause
        end
        clear thetablindstruct dist Zstruct
    end;
    
    if DISPLAY
        figure(17)
        %plot(phi,Thetablind,'+-b');
        %hold on
        plot(phi,Thetablindstruct,'-k');
        xlabel('Azimuth angle [deg]', labelOpts{:})
        ylabel('Horizon direction [deg]', labelOpts{:})
        pause
        close all
    end
end;
        
fov.BlindAngle = AntennaBlind;
fov.Phi = phi;

save(fovname,'fov')        
            
        
        
        
        




