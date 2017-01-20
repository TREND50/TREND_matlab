function [  ] = SelectRaw( nrun )
% Select events from list and write corresponding data to new file
% OMH 14/10/2013
% taken from SelectData.m

SharedGlobals

if SelectedData == 1
    disp 'Error! SharedGlobals::SelectedData == 1. Reading selected data instead of raw data! Fix this.'
    return
end
txtname = [TEXT_PATH sprintf('R%d_selection.txt',nrun)];
sel = load(txtname);


dets = sel(:,1);
evmin = sel(:,2);
evmax = sel(:,3);
nev = evmax-evmin+1;

for i = 1:length(dets)
        fclose all;
        % Open datafile for reading
        ft = OpenFileTime( nrun, dets(i) );
        fd = OpenFileData( nrun, dets(i) );
        if fd<0 | ft<0
            disp(sprintf('Could not find data for detector %d.',dets(i)))
            continue
        else
            disp(sprintf('Data file found for detector %d.',dets(i)))
        end

        % Open datafile for writing
        filename = sprintf( 'S%06d_A%04d_time.bin', nrun, dets(i) ); 
        ftw = fopen( [ SELDATA_PATH, filename ], 'w+' );
        filename = sprintf( 'S%06d_A%04d_data.bin', nrun, dets(i) ); 
        fdw = fopen( [ SELDATA_PATH, filename ], 'w+' );

        fseek( ft, 4*ibufft*(evmin(i)-1),'bof');
        fseek( fd, ibuff*(evmin(i)-1),'bof');
    
        TimeEvt = double( fread( ft,nev(i)*ibufft,'*uint32' ) );
        DataEvt = double( fread( fd, nev(i)*ibuff, '*uint8' ) );
    
        fwrite(ftw,TimeEvt,'*uint32');
        fwrite(fdw,DataEvt,'*uint8');

end
fclose all;
