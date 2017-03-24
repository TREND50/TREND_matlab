global C0 RAWDATA_PATH DST_PATH REFWE REFSN REFALT RAD2DEG DEG2RAD cutsettings labelOpts;

%addpath('../matlab_tools/')

CC=0;

%% Data Challenge or not
DC=0;
if DC == 1
  energy_eV = '5e17';
end
%% General
C0 = 2.99792458e8;  %m/s
DEG2RAD=pi/180;
RAD2DEG=180./pi;
phigeo = 182.66;  %Not to be trusted (referential?)
thetageo = 27.05;  %Not to be trusted (referential?)
% In GRAND ref (x=N, y=W, z=Up), direction pointed to:
phigeom = 2.66;  % West of North (<=>Declination = -2.66°)
thetageom = 152.92;  % Here theta = 0° <=> up. 152.92° <=> 62.95° below horizon (<=>Inclination = 62.95°)  
labelOpts = { 'FontSize', 14 };
scrsz = get( 0, 'ScreenSize' );

%% Site & DAQ
REFWE=0;
REFSN=0;
REFALT=2632;
PROPAGT = 0.77; %samples/m (=3.85ns/m)
FSAMPLING = 200e6;  %Hz
NBITS = 8;
SCALE = 3.3/2^NBITS;
NFFT = 256;
FREQMAX = 100e6;
FREQMIN = 50e6;
ErrorTrig = 2; % Trigger tagging error in units of samplee (<=> 10ns)
ErrorTrigScint = 4; % Trigger tagging error in units of samples (<=> 10ns)
ErrorAmp = 0.15; % Amplitue measurment error (fraction of total) 
runs2010 = [1:2400, 17886, 17959, 180234, 209212, 20979]; 
ALLDETS = [101:140 148:158];
ibuff = 1024;
ibuffs = 1024;
ibufft = 4;

%% Analysis

%% Preprocessing parameters

% Cuts for Candidates selection
cutsettings.T1.LimTheta = 85;
cutsettings.T2.NBoxesCut = 2;
cutsettings.T2.DurationCut = 2.5e-7;
cutsettings.T2.MinBoxes = 2;
cutsettings.T3.RatioMaxMin = 1.5;
%cutsettings.T4.LimDistBary = 300;  %Polar EW
cutsettings.T4.LimDistBary = 500;
cutsettings.T5.DistCut = 500;
HARD = 0;
if HARD == 1
    cutsettings.T5.time_span_T5_mn = 5; 
    cutsettings.T5.PhiCut = 10;
    cutsettings.T2.MinLong = 1;
    cutsettings.T5.NbComAntT5 = 3;
    cutsettings.T6.time_span_T6 = 60;
    cutsettings.T6.NbComAntT6 = 3;
    cutsettings.T6.NbEventsT6 = 2;
else
    cutsettings.T5.time_span_T5_mn = 3; %Polar NS
    cutsettings.T5.PhiCut = 5; %Polar NS
    cutsettings.T2.MinLong = 2;
    cutsettings.T5.NbComAntT5 = 4;
    cutsettings.T6.time_span_T6 = 30;
    cutsettings.T6.NbComAntT6 = 4;
    cutsettings.T6.NbEventsT6 = 3;
end
LIGHT = 0;  % Use light DSTs (no Candidate info)

% Elimination of events in a high trigger rate period
TrigRateFilterFlag=0;

% Elimination of consecutive coinc (fast) or not (slow)
QuickReject=1;

% Analysis type :
% 0: radio analysis
% 1: hybrid analysis
% 2: scintillator analysis
AnalysisType=0;

%%%%%%%%%%%% ONLY TRIGGER TIME BOX MODE ??? %%%%%%%%%%%%
%%%%%%%% only for antennas, impossible for scintillators
TrigTimeBox=1;

% Rejection of satured events in final dst building
FiltSat=0;

% Enable rejection through RawFilter
RawFilterFlag=1;

% Total Rejection of consecutive coincidences
ConsCoincCompleteRejection=1;

% Minimal multiplicity for coincidences (radio analysis)
RadioMiniMult=4;

% Record shower candidate signals
RecordSignals =0;

% Use selected data (Olivier's stuff)
SelectedData = 0;

% Remove empty signal soft (=0) or hard (=1) way
EmptySignalTotalFilt=1;

% Enable candidate selection
CandidateSelection=1;

% Use normal delay time or estimated delay time (for candidate selection)
UseDelay=0;

% Fast analysis: use time files only
AnaFast = 0;

%% Variables
CORREL = 1; % Perform cross correlation treatment
if AnalysisType>0
    CORREL = 0;
end
if AnaFast == 1
    CORREL = 0;
end
CORDELAYS = 0;
DISPLAY = 0;
TriggerTimeSpan = 60; %duration for trigger rate calculation [s]
TrigRateLimit = 100; % Cut periods with larger trigger rates [Hz]
msize=100000; % maximum size for structure

% RawFilter parameters, used in CoincidenceFiltering
threshold       = 5.0;
granularity_sampling = 15e-9; % s (!!!! used in RawFilter too !!!)
granularity_box  = 40e-9; % s (!!!! used in RawFilter too !!!)
max_out = 0;
max_out_ToT   = 300e-9; % s
max_block_ToT   = 350e-9; % s  % Was 350ns. Modified 28/02 for DAQproto events (large bounce)
noise_mode      = 1;
max_pretrig_ToT = 200e-9;

%% Background runs variables
size_event=1024; % bytes
nbevent=1024; % 1024 events of 1024 bytes by loop in background run

%% PATHS
if CC==0  % Local PCS, no environment variables defined
        RAWDATA_PATH = '../data/raw/';
        TEXT_PATH = './';
        CXX_PATH = '/home/martineau/TREND/cxx/';
        LOG_PATH = '../data/log/';
        MONITOR_PATH = '../data/monitor/';
        if DC
            CAND_PATH = '../data/candidates/candidates_dc/';
        else
            CAND_PATH = '../data/candidates/sauv/candidates_dst102014/';
        end
        PSD_PATH = '../data/psd/';
        STD_PATH = '../data/dst/dst_std/dst102014/';
        HYB_PATH = '../data/dst/dst_hyb/';
        SCI_PATH = '../data/dst/dst_scint/';
        if DC 
            DST_PATH = ['../data/dst/dst_dc/' energy_eV '/'];
        else
            if AnalysisType == 0
                DST_PATH = STD_PATH;
            elseif AnalysisType == 1
                DST_PATH = HYB_PATH;
            elseif AnalysisType == 2
                DST_PATH = SCI_PATH;
            end
        end
        BCKGR_PATH = '../data/back/';
        BACKDST_PATH = '../data/back/';
        SELDATA_PATH = '../data/raw/';
        CAL_PATH = '../data/cal/';
        SURVEY_PATH = RAWDATA_PATH;
        %SIMU_PATH = 'G:/data/simu/trend-50/1e17/';
        if HARD
            ACC_PATH = '../data/acceptance/hard/';
        else
            ACC_PATH = '../data/acceptance/soft/';
        end        
else  % @CC
    MONITOR_PATH = './';
    if DC  % Data challenge
        RAWDATA_PATH = [getenv( 'TREND_DATADC_PATH' )  energy_eV  '/'];
        if isempty( RAWDATA_PATH )
            disp( '$TREND_DATADC_PATH is undefined.' );
            RAWDATA_PATH = ['/sps/hep/trend/datachallenge/' energy_eV  '/'];
        end
    else  % Regular data
        RAWDATA_PATH = getenv( 'TREND_DATA_PATH' ); 
        if isempty( RAWDATA_PATH )
           disp( '$TREND_DATA_PATH is undefined.' );
  	       RAWDATA_PATH = '/sps/hep/trend/data/';
	    end
    end
    CXX_PATH = getenv( 'TREND_CXX_PATH' );
    if isempty( CXX_PATH )
        disp( '$TREND_CXX_PATH is undefined.' );
 	    CXX_PATH = '/afs/in2p3.fr/throng/trend/soft/ana/TREND_recons/';
    end
    TEXT_PATH = getenv( 'TREND_TEXT_PATH' );
    if isempty( TEXT_PATH )
        disp( '$TREND_TEXT_PATH is undefined.' );
        TEXT_PATH = '/afs/in2p3.fr/throng/trend/scratch/txt/';
    end
    if DC   %Data challenge
        DST_PATH = [getenv( 'TREND_DSTDC_PATH' ) energy_eV '/'];
        if isempty(DST_PATH)
            disp( '$TREND_DSTDC_PATH is undefined.' );
	        DST_PATH = ['/sps/hep/trend/dst_datachallenge/' energy_eV '/'];
        end
    else  % Regular data
        DST_PATH = getenv( 'TREND_DST_PATH' );       
        if isempty(DST_PATH)
            disp( '$TREND_DST_PATH is undefined.' );
            DST_PATH='/afs/in2p3.fr/home/throng/trend/scratch/dst/';
        end
    end
   if DC
       CAND_PATH = '/sps/hep/trend/TREND_candidates/dc/'
   else
       CAND_PATH =  '/sps/hep/trend/TREND_candidates/data/'
   end
    
    PSD_PATH = getenv( 'TREND_PSD_PATH' );
    if isempty( PSD_PATH )
        disp( '$TREND_PSD_PATH is undefined.' );
        PSD_PATH =   '/sps/hep/trend/psd/';
    end

    CAL_PATH = getenv( 'TREND_CAL_PATH' );
    if isempty( CAL_PATH )
        disp( '$TREND_CAL_PATH is undefined.' );
        CAL_PATH =   '/sps/hep/trend/dstcalib/';
    end
    %SIMU_PATH = '/sps/hep/trend/trend-50/';
    %LOG_PATH = '/sps/hep/trend/log/';
    %SELDATA_PATH = '/sps/trend/';
    %BACKDST_PATH = DST_PATH;
    %BCKGR_PATH = '/sps/hep/trend/data_bckgd/';
end

%% Filenames
coinc_filename='CoincData%d_%d.mat';
preprocess_filename='PreProcessData%d_%d.mat';
dst_filename='dst%d_%d.mat';
dsthyb_filename='dst%d_%d.mat'; % temporary
dstscint_filename='dstscint%d_%d.mat'; % temporary
calib_filename='calib%d_%d.mat';

%% Cut definition for EAS candidates
nAntsCut = 3; % no cut on antenna multiplicity
MinDistCut = 500; %in meters
DeltaCut = 10; % in degrees
ThetaCut = 90;
Chi2Cut = 10;
SlopeCut = 0.1;

% BarycenterRejection
DistBarySigma = 2.5;
DistBaryLimit=300;

% TrackFinder
TrackChi2PreCut = 50;
TrackTpsLimit = 61; % in seconds
TrackThetaLimit = 20; % in degrees
TrackPhiLimit = 20; % in degrees

