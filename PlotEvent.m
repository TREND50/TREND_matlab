function []=PlotEvent(Run,Ncoinc)


sharedGlobals;

evtsize=1024;
sampling=5e-9;
timescale=-((evtsize/2)-1)*sampling:sampling:(evtsize/2)*sampling;

%filename=sprintf('dst%d_%d.mat',Run,1);
filename=sprintf('PreProcessData%d_%d.mat',Run,1);
filename=[DST_PATH filename];
load(filename);


ind = find(Struct.Coinc.IdCoinc==Ncoinc)

Tag=Struct.Coinc.Det.Tag(ind,:);
%Time=Struct.Coinc.Det.Time(Ncoinc,:);
%TrigTime=[Struct.Coinc.Det.TrigTime(Ncoinc,:)];
Evt=Struct.Coinc.Det.Evt(ind,:);

Name=[Struct.Setup.Det.Name];
indsci=find([Struct.Setup.Det.isScint]==1);
indant=find([Struct.Setup.Det.isScint]==0);

indscitag=indsci(Tag(indsci)==1)
indanttag=indant(Tag(indant)==1)

% Caracs subplot
NbPlotFig=3;


NbPlotCurrent=1;
NbFig=1;

% Display options
DisplayBoxAnt=1;
DisplayBoxSci=1;



% Plot Antennas
if ~isempty(indanttag)
    
    dat=figure(1);
    scrsz=get(0,'ScreenSize');
    set(dat,'Position',[scrsz(1) scrsz(2) scrsz(3)*0.99 scrsz(4)*0.85]);
    
    for i=1:length(indanttag)
    
        % if more than plots than NbPlotFig, open a new fig
        if NbPlotCurrent>NbPlotFig
            NbPlotCurrent=1;
            NbFig=NbFig+1;
            dat=figure(NbFig);
            scrsz=get(0,'ScreenSize');
            set(dat,'Position',[scrsz(1) scrsz(2) scrsz(3)*0.99 scrsz(4)*0.85]);
        end;
    
        % Extraction data antenna
        ft = OpenFileData(Run,Name(indanttag(i)));
        if ft==-1
            display('Antenna not found');
            continue;
        end;
        start_event=(Evt(indanttag(i))-1)*evtsize+1;
        fseek(ft,start_event,'bof');
        data=fread(ft,evtsize,'uint8');
        fclose(ft);
    
        % subplotting event
        subplot(NbPlotFig,1,NbPlotCurrent)
        plot(timescale,data)
        eventtitle=sprintf('Antenna: %d ; Run : %d ; Nb Coinc : %d ; Ant Event : %d',Name(indanttag(i)),Run,Ncoinc,Evt(indanttag(i)));
        title(eventtitle);
        xlabel('Time (s)');
        ylim([0.9*min(data) 1.1*max(data)]);
        ylabel('ADC count');
        
        % Display signal boxes (coming from RawFilter)
        if DisplayBoxAnt
            
            hold on
            
            Sigma=Struct.Coinc.Det.Sigma(ind,indanttag(i));
            Mu=Struct.Coinc.Det.Mu(ind,indanttag(i));
            
            plot([timescale(1) timescale(end)],[Mu+4*Sigma Mu+4*Sigma],'k')
            plot([timescale(1) timescale(end)],[Mu-4*Sigma Mu-4*Sigma],'k')
            
            Carac=Struct.Coinc.Reject.CaracEvt(ind,indanttag(i),:);
            nbox=cell2mat(Carac(2));
            box_start=cell2mat(Carac(4));
            box_end=cell2mat(Carac(5));
            
            for k=1:nbox
                plot([timescale(1)+box_start(k)*sampling timescale(1)+box_start(k)*sampling],[0.95*min(data) 1.05*max(data)],'r')
                plot([timescale(1)+box_end(k)*sampling timescale(1)+box_end(k)*sampling],[0.95*min(data) 1.05*max(data)],'r')
                plot([timescale(1)+box_start(k)*sampling timescale(1)+box_end(k)*sampling],[0.95*min(data) 0.95*min(data)],'r')
                plot([timescale(1)+box_start(k)*sampling timescale(1)+box_end(k)*sampling],[1.05*max(data) 1.05*max(data)],'r')
            end;
            hold off
            
            clear CaracEvt nbox box_start box_end Sigma Mu
            
        end;
        
               
        clear data start_event title
    
        NbPlotCurrent=NbPlotCurrent+1;
    
    end;
end;

% Plot scintillators
if ~isempty(indscitag)

    NbFig=NbFig+1;
    dat=figure(NbFig);
    scrsz=get(0,'ScreenSize');
    set(dat,'Position',[scrsz(1) scrsz(2) scrsz(3)*0.99 scrsz(4)*0.85]);

    NbPlotCurrent=1;
    
    for i =1:length(indscitag)
        
        % Extraction data antenna
        ft = OpenFileData(Run,Name(indscitag(i)));
        if ft==-1
            display('Scintillator not found');
            continue;
        end;
        start_event=(Evt(indscitag(i))-1)*evtsize+1;
        fseek(ft,start_event,'bof');
        data=fread(ft,evtsize,'uint8');
        fclose(ft);
    
        % subplotting event
        subplot(NbPlotFig,1,NbPlotCurrent)
        plot(timescale,data)
        eventtitle=sprintf('Scintillator: %d ; Run : %d ; Nb Coinc : %d ; Scint Event : %d',Name(indscitag(i)),Run,Ncoinc,Evt(indscitag(i)));
        title(eventtitle);
        xlabel('Time (s)');
        ylim([0.8*min(data) 1.2*max(data)]);
        ylabel('ADC count');
        
        
        % Display signal boxes (coming from RawFilter)
        if DisplayBoxSci
            
            hold on
            
            Sigma=Struct.Coinc.Det.Sigma(ind,indscitag(i));
            Mu=Struct.Coinc.Det.Mu(ind,indscitag(i));
            
            plot([timescale(1) timescale(end)],[Mu+15*4*Sigma Mu+15*4*Sigma],'k')
            plot([timescale(1) timescale(end)],[Mu-15*4*Sigma Mu-15*4*Sigma],'k')
            
            Carac=Struct.Coinc.Reject.CaracEvt(ind,indscitag(i),:);
            nbox=cell2mat(Carac(2));
            box_start=cell2mat(Carac(4));
            box_end=cell2mat(Carac(5));
            
            nbox
            
            for k=1:nbox
                plot([timescale(box_start(k)) timescale(box_start(k))],[0.95*min(data) 1.05*max(data)],'r')
                plot([timescale(box_end(k)) timescale(box_end(k))],[0.95*min(data) 1.05*max(data)],'r')
                plot([timescale(box_start(k)) timescale(box_end(k))],[0.95*min(data) 0.95*min(data)],'r')
                plot([timescale(box_start(k)) timescale(box_end(k))],[1.05*max(data) 1.05*max(data)],'r')
            end;
            hold off
            
            clear CaracEvt nbox box_start box_end Sigma Mu
            
        end;
        clear data start_event title
    
        NbPlotCurrent=NbPlotCurrent+1;
        
    end;
    
end;
    