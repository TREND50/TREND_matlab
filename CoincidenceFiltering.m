function [NbIterDst]=CoincidenceFiltering(nrun,NbIter)


SharedGlobals;

NbCoinc=1;
NbIterDst=1;

RunTotalCoinc=msize;
TotalDet=80;

% Initialisation of final structure
FinalMult=zeros(RunTotalCoinc,1);
FinalMultAnt=zeros(RunTotalCoinc,1);
FinalMultSci=zeros(RunTotalCoinc,1);
FinalTag=zeros(RunTotalCoinc,TotalDet);
FinalMaxRaw=zeros(RunTotalCoinc,TotalDet);
FinalMinRaw=zeros(RunTotalCoinc,TotalDet);
FinalSat=zeros(RunTotalCoinc,TotalDet);
FinalStatus=-ones(RunTotalCoinc,TotalDet);
FinalFiltResult=zeros(RunTotalCoinc,TotalDet);
FinalIdCoinc=zeros(RunTotalCoinc,1);
FinalId=zeros(RunTotalCoinc,TotalDet);
FinalEvt=zeros(RunTotalCoinc,TotalDet);
FinalTime=zeros(RunTotalCoinc,TotalDet);
FinalUnixTime=zeros(RunTotalCoinc,TotalDet);
FinalTriggerRate=zeros(RunTotalCoinc,TotalDet);
FinalSigma=zeros(RunTotalCoinc,TotalDet);
FinalMu=zeros(RunTotalCoinc,TotalDet);

% Reloading RawFilter parameters
settings.granularity = granularity_box; % s
settings.max_out = max_out;
settings.max_out_ToT   = max_out_ToT; % s
settings.max_block_ToT   = max_block_ToT; % s
%settings.bounce_ratio    = bounce_ratio;
%settings.inhibit_window  = inhibit_window; % s
settings.noise_mode      = noise_mode;
settings.max_pretrig_ToT = max_pretrig_ToT;
%settings.pre_trig_win = pre_trig_win;

sumcoincraw=0;
sumcoincfilter=0;
sumok = 0;
sumtot = 0;

% Loop on sub dst
for i=1:NbIter
        
    % load sub structures
    filename=[DST_PATH sprintf(preprocess_filename,nrun,i)];
    load(filename);
     
    % Rename structure (Struct will be use for dst)
    PreStruct=Struct;
    clear Struct
    ncoinc=PreStruct.Setup.TotalCoinc;
    if ncoinc==0
       display(sprintf('No coinc found (Preprocess data in dst building stage) for run %d iteration %d',nrun,i));
        %NbDet=size(Tag,2);
	NbDet=length([PreStruct.Setup.Det.Name]);
	continue
    end;
    
    Tag=PreStruct.Coinc.Det.Tag;
    NbDet=size(Tag,2);
    
    Status = [PreStruct.Coinc.Det.Status];
    FiltResult = zeros(size(Status,1),size(Status,2));
    indsci=find([PreStruct.Setup.Det.isScint]==1);
    indant=find([PreStruct.Setup.Det.isScint]==0);
    Detectors = [PreStruct.Setup.Det.Name];
    %ncoinc=PreStruct.Setup.TotalCoinc;
    Carac=PreStruct.Coinc.Reject.CaracEvt;
    Sat=PreStruct.Coinc.Det.Sat;
    CoincId = PreStruct.Coinc.IdCoinc;
    Sigma = PreStruct.Coinc.Det.Sigma;
    Evt = PreStruct.Coinc.Det.Evt;
    sumcoincraw=sumcoincraw+ncoinc;
    
 
    %% COINCIDENCE REJECTION WITH RAWFILTER PARAMETERS
    
    % loop on coincidences   
    for j=1:ncoinc

        indscitag=indsci(Tag(j,indsci)==1);
        indanttag=indant(Tag(j,indant)==1);
        
        % Loop on antennas
        tagfilter=zeros(NbDet,1);
        tagfilter(indscitag)=1; % keep scintillators in any case
        %disp(sprintf('Checking signals for coinc %d',CoincId(j)))
        for k=1:length(indanttag)
            %disp(sprintf('Processing signal for antenna %d',Detectors(indanttag(k))))
            % RejectionRawFilter returns 1 if signal is OK, 0 if rejected.
            settings.threshold = (threshold+0)*Sigma(j,indanttag(k));
            if AnalysisType==0 && RawFilterFlag==1
                [tagfilter(indanttag(k)) filtresult]=RejectionRawFilter(Carac(j,indanttag(k),:),settings);
                
            else
                tagfilter(indanttag(k))=1; %keep the value is RawFilter is disabled (hybrid analysis)
            end;
            
            if tagfilter(indanttag(k))== 0 
                Status(j,indanttag(k)) = Status(j,indanttag(k))+2;  % Bit1 set at 1
                ind=find(filtresult(:,2)==1);
                FiltResult(j,indanttag(k))=ind;
            end
            
            if FiltSat==1 & Sat(j,indanttag(k))==1
                tagfilter(indanttag(k))=0;
            end;
        end;
        
        tagfinal= tagfilter;
        multfilter=length(find(tagfinal==1));
		multsci=length(find(tagfinal(indsci)==1));
        if multfilter>3 || multsci==3
            
            indtagfinal=find(tagfinal==1);  
            FinalMult(NbCoinc)=multfilter;
            FinalMultAnt(NbCoinc)=length(find(tagfinal(indant)==1));
            FinalMultSci(NbCoinc)=length(find(tagfinal(indsci)==1));
            FinalTag(NbCoinc,1:NbDet)=tagfinal;
            FinalMinRaw(NbCoinc,1:NbDet)=PreStruct.Coinc.Det.MinRaw(j,1:NbDet);
            FinalMaxRaw(NbCoinc,1:NbDet)=PreStruct.Coinc.Det.MaxRaw(j,1:NbDet);
            FinalSat(NbCoinc,1:NbDet)=Sat(j,1:NbDet);
            FinalStatus(NbCoinc,1:NbDet)=Status(j,1:NbDet);
            FinalFiltResult(NbCoinc,1:NbDet)=FiltResult(j,1:NbDet);
            FinalIdCoinc(NbCoinc)=PreStruct.Coinc.IdCoinc(j);
            FinalId(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.Id(j,indtagfinal);
            FinalEvt(NbCoinc,1:NbDet)=PreStruct.Coinc.Det.Evt(j,1:NbDet);
            FinalTime(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.Time(j,indtagfinal);
            FinalUnixTime(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.UnixTime(j,indtagfinal);
            FinalTriggerRate(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.TriggerRate(j,indtagfinal);
            FinalSigma(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.Sigma(j,indtagfinal);
            FinalMu(NbCoinc,indtagfinal)=PreStruct.Coinc.Det.Mu(j,indtagfinal);
            
            NbCoinc=NbCoinc+1;
            
            if NbCoinc>msize  % Write to dst now
                
                NbCoinc=NbCoinc-1;
                sumcoincfilter=sumcoincfilter+NbCoinc;

                Setup=PreStruct.Setup;
                Struct.Setup=Setup;
                Struct.Setup.TotalCoinc=NbCoinc;
                Struct.Coinc.Mult=FinalMult(1:NbCoinc);
                Struct.Coinc.MultAnt=FinalMultAnt(1:NbCoinc);
                Struct.Coinc.MultSci=FinalMultSci(1:NbCoinc);
                Struct.Coinc.IdCoinc=FinalIdCoinc(1:NbCoinc);
                Struct.Coinc.Det.Tag=FinalTag(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.MaxRaw=FinalMaxRaw(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.MinRaw=FinalMinRaw(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Sat=FinalSat(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Status=FinalStatus(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.FiltResult=FinalFiltResult(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Id=FinalId(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Evt=FinalEvt(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Time=FinalTime(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.UnixTime=FinalUnixTime(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.TriggerRate=FinalTriggerRate(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Sigma=FinalSigma(1:NbCoinc,1:NbDet);
                Struct.Coinc.Det.Mu=FinalMu(1:NbCoinc,1:NbDet);


                filename=[DST_PATH sprintf(dst_filename,nrun,NbIterDst)];
                save(filename,'Struct');
                clear Struct
                NbIterDst=NbIterDst+1;
                NbCoinc=1;
                
                % Initialisation of final structure
                FinalMult=zeros(RunTotalCoinc,1);
                FinalMultAnt=zeros(RunTotalCoinc,1);
                FinalMultSci=zeros(RunTotalCoinc,1);
                FinalTag=zeros(RunTotalCoinc,TotalDet);
                FinalMaxRaw=zeros(RunTotalCoinc,TotalDet);
                FinalMinRaw=zeros(RunTotalCoinc,TotalDet);
                FinalSat=zeros(RunTotalCoinc,TotalDet);
                FinalStatus=-ones(RunTotalCoinc,TotalDet);
                FinalFiltResult=zeros(RunTotalCoinc,TotalDet);
                FinalIdCoinc=zeros(RunTotalCoinc,1);
                FinalId=zeros(RunTotalCoinc,TotalDet);
                FinalEvt=zeros(RunTotalCoinc,TotalDet);
                FinalTime=zeros(RunTotalCoinc,TotalDet);
                FinalUnixTime=zeros(RunTotalCoinc,TotalDet);
                FinalTriggerRate=zeros(RunTotalCoinc,TotalDet);
                FinalSigma=zeros(RunTotalCoinc,TotalDet);
                FinalMu=zeros(RunTotalCoinc,TotalDet);
                
            end;          
            
            clear indtagfinal
            
        end;
        clear mult indscitag indanttag tagfilter multfilter
        
        decim = floor(ncoinc/10);
        if floor(j/decim)==j/decim
            display(sprintf('Signal filtering done at %2.0f percent for PreProcessData_%d ',j/decim*10,i));
        end        
    end;

    for k=1:length(Detectors)
      indc = find(Tag(:,k)==1);
      %CoincId(indc)
      res = FiltResult(indc,k);
      nok = length(find(res==0));
      nout = length(find(res==1));
      npre = length(find(res==2));
      ncen = length(find(res==3));      
      nrep = length(find(res==4));
      nnob = length(find(res==5));    
      nno = nout+npre+ncen+nrep+nnob;
      ntot = length(indc);
      disp(sprintf('Antenna %d in %d coincs: %d events OK (%3.1f pc), %d rejected (%3.1f pc)',Detectors(k),ntot,nok,nok/ntot*100,nno,nno/ntot*100));
      disp(sprintf('%d events rejected: %d too long central; %d too early central; %d too long outside, %d pulses outside, %d no box',nno,ncen,npre,nout,nrep,nnob))
      % dump list to file
      filename = [TEXT_PATH sprintf('R%d_A%d_filtering.mat',nrun,Detectors(k))];
      good = Evt(indc(res==0),k);
      bad = Evt(indc(res>0),k);
      save(filename,'good','bad')
      sumok = sumok+nok;
      sumtot = sumtot+ntot;
      %pause
    end
end

NbCoinc=NbCoinc-1;
sumcoincfilter=sumcoincfilter+NbCoinc;
Setup=PreStruct.Setup;
Struct.Setup=Setup;
Struct.Setup.TotalCoinc=NbCoinc;
Struct.Coinc.Mult=FinalMult(1:NbCoinc);
Struct.Coinc.MultAnt=FinalMultAnt(1:NbCoinc);
Struct.Coinc.MultSci=FinalMultSci(1:NbCoinc);
Struct.Coinc.IdCoinc=FinalIdCoinc(1:NbCoinc);
Struct.Coinc.Det.Tag=FinalTag(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.MaxRaw=FinalMaxRaw(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.MinRaw=FinalMinRaw(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Sat=FinalSat(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Status=FinalStatus(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.FiltResult=FinalFiltResult(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Id=FinalId(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Evt=FinalEvt(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Time=FinalTime(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.UnixTime=FinalUnixTime(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.TriggerRate=FinalTriggerRate(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Sigma=FinalSigma(1:NbCoinc,1:NbDet);
Struct.Coinc.Det.Mu=FinalMu(1:NbCoinc,1:NbDet);

filename=[DST_PATH sprintf(dst_filename,nrun,NbIterDst)];
save(filename,'Struct');
                
display(sprintf('Total coinc found: %d',sumcoincraw))
display(sprintf('Total coinc remaining: %d',sumcoincfilter))
display(sprintf('Filter efficiency (coincs): %3.1f pc',100*(1-(sumcoincfilter/sumcoincraw))))
display(sprintf('Filter efficiency (signals): %3.1f pc',100*(1-(sumok/sumtot))))



  
