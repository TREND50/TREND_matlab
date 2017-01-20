function [ ft, filename ] = OpenFileBackgroundTime( nrun, antenna )
% Ouvre le fichier time correspondant � l'antenne rep�r�e par le
% num�ro antenna, dans le run rep�r� par le num�ro nrun.

SharedGlobals

BCKGR_PATH=RAWDATA_PATH;

filename = sprintf( 'R%06d_A%04d_BACK_time.bin', nrun, antenna ); 
ft = fopen( [ BCKGR_PATH, filename ], 'r' );
if ft == -1
  filename = sprintf( 'B%04d_%d_time.bin', nrun, antenna );
  ft = fopen( [ BCKGR_PATH, filename ], 'r' );
end
if ft ==-1
  filename = sprintf( 'B%04d_u%d_time.bin' ,nrun, antenna );
  ft = fopen( [ BCKGR_PATH, filename ], 'r' );
end
if ft ==-1
  filename = sprintf( 'B%05d_%d_time.bin' ,nrun, antenna );
  ft = fopen( [ BCKGR_PATH, filename ], 'r' );
end
