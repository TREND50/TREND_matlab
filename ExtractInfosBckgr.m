function [nb_loop,ant,big]=ExtractInfosBckgr(nrun)
% Extrait le nbre de loops et les antennes utilisés dans le run background

min_ant = 100; % valeur à modifier en fonction de la config
max_ant=158; % valeur à modifier en fonction de la config
ant=[];
nb_loop=[];
big=0;  % vérification de la taille du fichier

ft = OpenFileBackgroundTime(nrun,min_ant);
while (ft<0) & (min_ant<max_ant+1)
    min_ant = min_ant + 1;
    ft = OpenFileBackgroundTime(nrun,min_ant);
end
   
if min_ant <= max_ant % At least one file was found for this run

  fclose( ft );
  for nant = 1:( max_ant - min_ant + 1 );
    antenna = min_ant + nant - 1; % true antenna ID
    %fd = OpenFileBackground(nrun,antenna);
    ft = OpenFileBackgroundTime(nrun,antenna);
    if (ft ~= -1) 
      ant=[ant antenna];
      %data=fread(fd);    
      time=fread(ft,inf,'uint32');
      nbloop=length(time)/4;
      nb_loop=[nb_loop nbloop];
      %nbevent = [nbevent length(data)/(1024*nbloop)]; % 1024 bytes dans les runs background
      fclose(ft);
      %fclose (fd);
    end
  end
end
if max(nb_loop)>30
    big=1;
end;
