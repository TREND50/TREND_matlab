function [] = ReadFetchFile(nrun)
% Function to get waveforms from teh file fetchR[nrun].txt
% OMH 21/06/2013

SharedGlobals;

%% Read list of waveforms
listname = sprintf('~/fetchR%d.txt',nrun);
if fopen(listname)<0
    disp(sprintf('Could not find %s',listname))
    return
end
list = load(listname);
antid = list(:,1);
evtid = list(:,2);
coincid = list(:,3);

%% Loop on waveforms, write them to file.
for i=1:length(antid)
   disp(sprintf('Selecting event %d for antenna %d',evtid(i),antid(i))) 
   fd = OpenFileData( nrun, antid(i));
   if fd>0
       startpoint=(evtid(i)-1)*ibuff+1;
       fseek(fd,startpoint,'bof');
       DataEvt = double( fread( fd, ibuff, '*uint8' ) );
       filename = sprintf( 'R%06d_A%04d_C%d_data.bin', nrun,antid(i),coincid(i) ); 
       fdw = fopen(filename, 'w+' );
       fwrite(fdw,DataEvt,'*uint8');
       fclose(fdw)
   else
       disp(sprintf('Could not find data file for R%d A%d...',nrun,antid(i)))
   end
end
