function []=StructureBuilder(nrun)
  

nrun=sscanf(nrun,'%d')

% load global variables
SharedGlobals;

if DC == 1 % Data challenge
        display(sprintf('DataChallenge = %d, setting DST_PATH to %s, RAWDATA_PATH to %s',DC,DST_PATH,RAWDATA_PATH)) 
end

% Check for existing raw dst
% - if only CoincData, start from pre-analysis
% - if only PreProcessData, start from filtering
% - if dst, ask for rebuilding (from which step)
% - if nothing, start from beginning!

% Define Building flags
DstBuild=0;
PreProcessBuild=0;
CoincBuild=0;

NbIter=1; % just looking for the first file
filename = [DST_PATH sprintf(dst_filename,nrun,NbIter)];
ft=fopen(filename);
filename2 = [DST_PATH sprintf(preprocess_filename,nrun,NbIter)];
ft2=fopen(filename2);
filename3 = [DST_PATH sprintf(coinc_filename,nrun,NbIter)];
ft3=fopen(filename3);

if (ft~=-1) % dst file found
  
    display('Dst found for this run')
    display('Dst will be rebuilt from preprocessed data')
    display('Previous dst will be erased')
    fclose(ft);
    pause(5);
    DstBuild=1; % only dst rebuild
    
elseif ft2~=-1 % preprocess file found
    
    display('Preprocess data found for this run')
    display('Dst will be built from preprocessed data')
    fclose(ft2);
    DstBuild=1;
    
elseif ft3~=-1 % coinc file found
    
    display('Coinc data found for this run')
    display('Preprocess data will be generated')
    display('Dst will be built from preprocessed data')
    fclose(ft3);
    PreProcessBuild=1;
    DstBuild=1;
    
else
    % No file found, generate structure from beginning
    CoincBuild=1;
    PreProcessBuild=1;
    DstBuild=1;
    
end;

% if we don't start from the beginning, we need the number of iteration
if CoincBuild==0
    
    filename = [DST_PATH sprintf(coinc_filename,nrun,NbIter)];
    ft=fopen(filename);
    while ft~=-1
        NbIter=NbIter+1; % initialised at 1 previously
        filename = [DST_PATH sprintf(coinc_filename,nrun,NbIter)];
        ft=fopen(filename);
        fclose all;
    end;
    NbIter=NbIter-1; % final value
    
end;



%%%%%%%%%%%%%%%%%%% PRE PROCESSING DATA OPERATIONS%%%%%%%%%%%%%%%%%%%%%%%%


% Coincidence building
if CoincBuild
    
    display('Building setup configuration...')

    [Struct]=RunSetupBuilder(nrun,AnalysisType);
    RunSetup = Struct.Setup;
    if ~exist('RunSetup')
        display('Setup building error...')
        return
    else
        display('Done!')
    end

    if Struct.Setup.TotalEvt>100e6
        disp(sprintf('%d triggers in run, exceeds limit! Aborting.',Struct.Setup.TotalEvt));
		return
    else
        % OBSOLETE: replaced by (independant) CalibStructureBuilder
        % 18/05/2012 OMH
        %display('Loading calibration data...')
        %[Struct]=CalibBuilder(Struct);
        %display('Done!')

        display('Building Time Table...')
        [EventTimeTable,Struct]=EventTimeTableBuilder(Struct,TrigRateFilterFlag);
        display('Done!')
        save('/afs/in2p3.fr/throng/trend/soft/ana/TREND_insertsimu/table.mat','EventTimeTable')

        if (Struct.Setup.TotalEvt-size(Struct.Setup.InfosRun.EmptySignals,1))==0
           display(sprintf('All events empty for run %d',nrun));
            filename = sprintf(coinc_filename,nrun,1);
            EmptyStructure(filename,Struct,nrun,1);
            filename = sprintf(preprocess_filename,nrun,1);
            EmptyStructure(filename,Struct,nrun,1);
            filename = sprintf(dst_filename,nrun,1);
            EmptyStructure(filename,Struct,nrun,1);
            return
        end;
        
        % Saving total structure
        StructBase=Struct;

        % Creation of sub-structures of msize events
        NbIter=0;
        StopOrder=0;
        RefTable=1;

        while StopOrder==0

            NbIter=NbIter+1;       
            % Saving filename
            filename = [DST_PATH sprintf(coinc_filename,nrun,NbIter)];      
            % Reduction of EventTimeTable (depending on number of iterations)
            EventTimeTable=EventTimeTable(RefTable:end,:);       
            % Look for coincidences in the reduced time table, give the new RefTable
            [Struct,RefTable]=CoincidenceFinder(StructBase,EventTimeTable,NbIter,AnalysisType); 
            [GlobalCoincRate,DetCoincRate]=CoincidenceRate(Struct);
            Struct.Setup.InfosRun.GlobalCoincRateRaw=GlobalCoincRate;
            Struct.Setup.InfosRun.DetCoincRateRaw=DetCoincRate;
            % Saving substructures
            save(filename,'Struct');
            clear Struct;

            % If RefTable is null, there is no more event to study
            if RefTable==0
                StopOrder=1;
            end;
        end;		
    end;
end;


% Pre process building
if PreProcessBuild
    
    % Loop on all files
    for i=1:NbIter
    
        filename =[DST_PATH sprintf(coinc_filename,nrun,i)];
        load(filename);
        
        if Struct.Setup.TotalCoinc==0
			display(sprintf('No coinc found (CoincData stage) for run %d iteration %d',nrun,i));
            filename = sprintf(preprocess_filename,nrun,i);
            EmptyStructure(filename,Struct,nrun,i);
            display(sprintf('DST %s now saved to file.',filename));
            continue
        end

        % Consecutive coincs filter
        display('Lauching Consecutive coincidences filter...');
        [Struct]=ConsecutiveCoincidenceFilter(Struct);
        display('Done!')
        
        if QuickReject & SelectedData == 0 % structure modification: erase rejected events
            goodevts=find([Struct.Coinc.Reject.ConsCoinc]==0);
	   
            if ~isempty(goodevts)

                Struct.Setup.TotalCoinc=length(goodevts);
                Struct.Coinc.Mult=Struct.Coinc.Mult(goodevts);
                Struct.Coinc.MultAnt=Struct.Coinc.MultAnt(goodevts);
                Struct.Coinc.MultSci=Struct.Coinc.MultSci(goodevts);
                Struct.Coinc.IdCoinc=Struct.Coinc.IdCoinc(goodevts);
                Struct.Coinc.Det.Tag=Struct.Coinc.Det.Tag(goodevts,:);
                Struct.Coinc.Det.Status=Struct.Coinc.Det.Status(goodevts,:);
                Struct.Coinc.Det.Id=Struct.Coinc.Det.Id(goodevts,:);
                Struct.Coinc.Det.Evt=Struct.Coinc.Det.Evt(goodevts,:);
                Struct.Coinc.Det.Time=Struct.Coinc.Det.Time(goodevts,:);
                Struct.Coinc.Det.UnixTime=Struct.Coinc.Det.UnixTime(goodevts,:);
                Struct.Coinc.Det.TriggerRate=Struct.Coinc.Det.TriggerRate(goodevts,:);
            
                [GlobalCoincRate,DetCoincRate]=CoincidenceRate(Struct);
                Struct.Setup.InfosRun.GlobalCoincRateQuickReject=GlobalCoincRate;
                Struct.Setup.InfosRun.DetCoincRateQuickReject=DetCoincRate;
	
            else
                display(sprintf('No coinc found (Quick reject stage) for run %d iteration %d',nrun,i));
                filename = sprintf(preprocess_filename,nrun,i);
                EmptyStructure(filename,Struct,nrun,i);
                display(sprintf('DST %s now saved to file.',filename));
                continue
            end;            
        end;
        
        % Raw filter (Valentin)
        display('Lauching raw filter...');
        [Struct] = RawFilter(Struct);
        display('Done!')
        
        % Correction of scintillators trigger time using RawFilter result    
        % SciTimeCorrection2: use alternative method (strong deviation from mean)
        [Struct]=SciTimeCorrection2(Struct);
        
        % Saving pre process data
        filename = [DST_PATH sprintf(preprocess_filename,nrun,i)];
        save(filename,'Struct');
        clear Struct
    
    end
    
end;



%%%%%%%%%%%%%%%%%%%%%%%%% BUILDING DST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DstBuild
    
    display('Building final DST...')
    
    [NbIterDst]=CoincidenceFiltering(nrun,NbIter);
    
    for i=1:NbIterDst
        
        filename =[DST_PATH sprintf(dst_filename,nrun,i)];
        load(filename);
        
        if Struct.Setup.TotalCoinc==0
            display(sprintf('No coinc found (dst building stage) for run %d iteration %d',nrun,i));
            filename = sprintf(dst_filename,nrun,i);
            EmptyStructure(filename,Struct,nrun,i);
            continue
        end;
        
        display('Launching TriggerTimeBuilder...')
        [Struct]=TriggerTimeBuilder(Struct);
        display('Done!')
        
        display('Launchig StatusBuilder...')
        [Struct]=StatusBuilder(Struct);
        display('Done!')
        
        % Perform cross correlation treatment
        if CORREL
            display('Launching CrossCorrelationBuilder...')
            [Struct]=CrossCorrelationBuilder(Struct);
        else
            Struct.Coinc.Det.TrigCor = zeros(Struct.Setup.TotalCoinc,length([ Struct.Setup.Det.Name ]));
            Struct.Coinc.Det.CoefCor = zeros(Struct.Setup.TotalCoinc,length([ Struct.Setup.Det.Name ]));
        end
        
        % Build .txt file for C++ reconstruction (fast)
        display('Now reconstructing pure radio coincs only...')
        display('Launching CoinctableFileBuilder...')       
        [Struct]=CoinctableFileBuilder(Struct,0); 
        display('Done!')
        reconsex = [CXX_PATH 'recons'];
        unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d',TEXT_PATH,reconsex,nrun,CORREL));
        %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d',TEXT_PATH,reconsex,nrun,CORREL));
        [Struct]=ReconsLoader(Struct,0);
        [Struct]=Dist2Source(Struct); % distance antennas - point sources
        display('Launching CandidateFileBuilder...')
        [Struct]=CandidateFileBuilder(Struct);
        display('Done!')

        %reconsShowerex = [CXX_PATH 'reconsShower'];
        %%unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,0));
        %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL, 0));
        [Struct]=ShowerLoader(Struct,0);
        %%unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,1));
        %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,1));
        [Struct]=ShowerLoader(Struct,1);
        
        display('Now reconstructing hybrid coincs only...')
        display('Launching CoinctableFileBuilder...')       
        [Struct]=CoinctableFileBuilder(Struct,1); 
        display('Done!')
        %reconsScintex = [CXX_PATH 'reconsScint'];
        %unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d',TEXT_PATH,reconsScintex,nrun,CORREL));
        %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d',TEXT_PATH,reconsScintex,nrun,CORREL));
        [Struct]=ReconsLoader(Struct,1);        
        display('Done!')

        % Recalculation for corrected delays       
        if 0
	   Struct.Coinc.DelayCorrRecons=[];
 	   display('Recalculation for new delays')
 	   display('Now reconstructing pure radio coincs only...')
 	   display('Launching CoinctableFileBuilder...')
 	   [Struct]=CoinctableFileBuilder(Struct,0);
 	   display('Done!')
 	   reconsex = [CXX_PATH 'recons'];
 	   unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d',TEXT_PATH,reconsex,nrun,CORREL));
 	   %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d',TEXT_PATH,reconsex,nrun,CORREL));
 	   [Struct]=ReconsLoader(Struct,0);
 	   [Struct]=Dist2Source(Struct); % distance antennas - point sources
 	   display('Launching CandidateFileBuilder...')
 	   [Struct]=CandidateFileBuilder(Struct);
 	   display('Done!')
 	   reconsShowerex = [CXX_PATH 'reconsShower'];
 	   unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,0));
 	   %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL, 0));
 	   [Struct]=ShowerLoader(Struct,0);
 	   unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,1));
 	   %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d %d',TEXT_PATH,reconsShowerex,nrun,CORREL,1));
 	   [Struct]=ShowerLoader(Struct,1);
 
 	   display('Now reconstructing hybrid coincs only...')
 	   display('Launching CoinctableFileBuilder...')
 	   [Struct]=CoinctableFileBuilder(Struct,1);
 	   display('Done!')
 	   %reconsScintex = [CXX_PATH 'reconsScint'];
 	   %unix(sprintf('export TREND_TEXT_PATH=%s;%s %d %d',TEXT_PATH,reconsScintex,nrun,CORREL));
 	   %unix(sprintf('setenv TREND_TEXT_PATH %s;%s %d %d',TEXT_PATH,reconsScintex,nrun,CORREL));
 	   [Struct]=ReconsLoader(Struct,1);
 	   display('Done!')
        end
        
        % NO SCINT RECONS
        display('Now reconstructing pure scintillator coincs only...')
        [Struct]=SciOnlyRecons(Struct, AnalysisType);
        
        % Structure save       
        %Struct.Coinc.IsShower=4
        filename =[DST_PATH sprintf(dst_filename,nrun,i)];
        save(filename,'Struct');
        display(sprintf('DST %s now saved to file.',filename));
        clear Struct
        
    end;
end;

    
    
