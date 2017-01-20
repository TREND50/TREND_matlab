function [  ] = SelectData( nrun )
% Select events from list and write corresponding data to new file
% OMH 27/12/2011

SharedGlobals;

dstname = [DST_PATH sprintf('R%d_sel.mat',nrun)];
dst = load(dstname);


sel = sortrows(dst.selection,2);
icoinc = sel(:,1);
idet = sel(:,2);
ievt = sel(:,3);
ic = sort(icoinc);
ncoinc = sum(diff(ic)>0)
nevents = length(icoinc)
pause
id = 0;
for i = 1:nevents
    if id~=idet(i)  % new detector
        fclose all;
        % Open datafile for reading
        ft = OpenFileTime( nrun, idet(i) );
        fd = OpenFileData( nrun, idet(i) );
        if fd<0 | ft<0
            disp(sprintf('Could not find data for detector %d.',idet(i)))
            continue
        end
        id = idet(i);
        j = 1;
        % Open datafile for writing
        filename = sprintf( 'S%06d_A%04d_time.bin', nrun, id ); 
        ftw = fopen( [ SELDATA_PATH, filename ], 'w+' );
        filename = sprintf( 'S%06d_A%04d_data.bin', nrun, id ); 
        fdw = fopen( [ SELDATA_PATH, filename ], 'w+' );

    end
    fseek( ft, 4*ibufft*(ievt(i)-1),'bof');
    fseek( fd, ibuff*(ievt(i)-1),'bof');
    
    TimeEvt = double( fread( ft,ibufft,'*uint32' ) );
    DataEvt = double( fread( fd, ibuff, '*uint8' ) );
    
    fwrite(ftw,TimeEvt,'*uint32');
    fwrite(fdw,DataEvt,'*uint8');
    
    ipos(i) = j;
    j = j+1;

end
fclose all;

selection = [icoinc idet ievt ipos'];
save(dstname,'selection')