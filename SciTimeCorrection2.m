function [Struct]=SciTimeCorrection2(Struct)

 SharedGlobals; 

 RunSetup = Struct.Setup;
 nrun = RunSetup.Run;
 Ncoinc = RunSetup.TotalCoinc;
 indsci = find([RunSetup.Det.isScint]==1);
 Det = [RunSetup.Det.Name];
 Tag = [Struct.Coinc.Det.Tag];
 Time = [Struct.Coinc.Det.Time];
 Evt = [Struct.Coinc.Det.Evt];
 
 TimeSave = Time;
% DISPLAY = 1

 for i=1:length(indsci)
     fdc = OpenFileData( nrun, Det(indsci(i)));
     ftc = OpenFileTime( nrun, Det(indsci(i)));
     indtag = find(Tag(:,indsci(i))==1);
     ev = Evt(indtag,indsci(i));
     
     for j=1:length(indtag)
         
            fseek(fdc,(ev(j)-1)*ibuffs,'bof'); % Skip events
            buff=fread(fdc,ibuffs,'*uint8');
            
            fseek(ftc,(ev(j)-1)*ibufft*4,'bof'); % Skip events
            buff2=fread(ftc,ibufft,'*uint32');
            dbuff2=double(buff2);
            if ismember(nrun,runs2010) % Old soft
                trigpos = dbuff2(4);    
            else  % New soft
                trigpos = ibuffs/2;
            end
            
            v = double(buff);
            vfilt = v-mean(v);
            [vmax imax] = max(vfilt);
            if imax>=ibuffs/2  % Pulse in the end
                samp = [1:200];
            else
                samp = [ibuffs-200:ibuffs];
            end
            %samp = [1:200];
            stdclean = std(vfilt(samp));  % Std for 0 signal
            mv = mean(vfilt);
            cut = mv+10*stdclean;
            sig = find(abs(vfilt)>cut & abs(vfilt-mv)>10);
            if size(sig,1)>0
                trigposcor = sig(1);
            else
                trigposcor = trigpos;
            end
            delta(j) = trigposcor-trigpos;
            
            % Boxes from RawFilter
         %Carac=Struct.Coinc.Reject.CaracEvt(indtag(j),indsci(i),:);
         %nbox=cell2mat(Carac(2));
         %box_start=cell2mat(Carac(4));
         %box_amp=cell2mat(Carac(7));
         %
         %indamp=find(box_amp>4);
         %if ~isempty(indamp)
         %   shiftdeltat(j)=(box_start(indamp(1))-(ibuffs/2));
         %else
         %    shiftdeltat(j)=(box_start(1)-(ibuffs/2));
         %end;
         Time(indtag(j),indsci(i))=Time(indtag(j),indsci(i))+delta(j);
         %Time(indtag(j),indsci(i))=Time(indtag(j),indsci(i))+shiftdeltat(j);

         if DISPLAY 
            figure(1) 
            hold off
            ind =[1:ibuffs]; 
            plot(ind,v,'k')
            hold on
            plot(ind,abs(vfilt),'r')
            trigline = line([trigpos trigpos],[min(v)*0.5 max(v)*1.5]);
            set(trigline,'LineWidth',2,'Color','k');
            corline2 = line([trigposcor trigposcor],[0  255]);
            set(corline2,'LineWidth',2,'Color','r');
            trigposcor1 = ibuffs/2+shiftdeltat(j);
            corline1 = line([trigposcor1 trigposcor1],[0  255]);
            set(corline1,'LineWidth',2,'Color','b');
            cutline = line([1 ibuffs],[cut cut]);
            set(cutline,'LineWidth',1,'Color','g');
            xlabel(sprintf('Scintillator %d Evt %d',Det(indsci(i)),ev(j)),labelOpts{:})
            [maxi imaxi] = max(v);
            %xlim([imaxi-100 imaxi+200])
            ylim([0 255])
            bline = mean(v(1:500));
            integ = cumtrapz(abs(v-bline));
            plot(ind,integ,'b');

            hold off
         end
 
         if DISPLAY
            delta(j)
            shiftdeltat(j)
            pause
        end
         
     end;
 end 
 
 if DISPLAY
    difcor = delta-shiftdeltat;
    figure(2)
    hist(difcor,100)
 end
 
 Struct.Coinc.Det.Time=Time;
 Struct.Coinc.Det.TimeSave=TimeSave;

