function [ fd, filename ] = OpenFileBackground( nrun, antenna )
% Ouvre le fichier BCKGR correspondant a l'antenne reperee par le
% numero ant, dans le run repere par le numero nrun.

SharedGlobals
%BCKGR_PATH
%BCKGR_PATH=RAWDATA_PATH;

filename = sprintf( 'R%06d/R%06d_A%04d_BACK_data.bin', nrun,nrun, antenna );
fd = fopen( [ BCKGR_PATH, filename ], 'r' );
if fd == -1
    filename =sprintf( 'R%06d_A%04d_BACK_data.bin',nrun, antenna );
    fd = fopen( [ BCKGR_PATH, filename ], 'r' );
end
if fd == -1
    filename = sprintf( 'B%04d_%d_data.bin', nrun, antenna );
    fd = fopen( [ BCKGR_PATH, filename ], 'r' );
end
if fd == -1
  filename = sprintf( 'B%04d_u%d_data.bin', nrun, antenna );
  fd = fopen( [ BCKGR_PATH, filename ], 'r' );
end
if fd == -1
  filename = sprintf( 'B%05d_%d_data.bin', nrun, antenna );
  fd = fopen( [ BCKGR_PATH, filename ], 'r' );
end
