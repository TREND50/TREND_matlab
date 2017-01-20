function [] = loopDummyEvent()
% Run generateDummyEvent and plot resulting resiolution
% OMH 15/02/2014


SharedGlobals;

%% Init
ntries = 1000
rhoi = 1e8;  % Distant source: plane = sph recons
phii = zeros(1,ntries);
thetai = zeros(1,ntries);

%% Loop
for i = 1:ntries
    phii(i) = 360*rand(1,1);
    cos_th = rand(1,1);   
    thetai(i) = acosd(cos_th);
    generateDummyEvent(rhoi,thetai(i),phii(i),i)
    fclose all;
end

disp 'Now perform C++ recons...'
pause

for i = 1:ntries
  %% Load files
  filename = [TEXT_PATH sprintf( 'R%d_planerecons.txt',i)];
  if fopen( filename )<0
    disp 'No recons file. Abort.'
    return
  end
  planRes = load( filename );
  %
  filename = [TEXT_PATH sprintf('R%d_sphrecons.txt', i)];
  if fopen(filename)<0
      disp(sprintf('File %s does not exist.',filename));
      return
  end
  sphRes = load( filename );
  %% Load variables
  timep = planRes(1,2);
  multp = planRes(1,3);
  thetap(i) = planRes(1,4);
  dThetap(i) = planRes(1,5);
  phip(i) = planRes(1,6);
  dPhip(i) = planRes(1,7);
  chi2p(i) = planRes(1,8);
  signif = planRes(1,9);
  thetap(i) = 180-thetap(i);
  phip(i) = phip(i)+180;
  phip(i)= mod(phip(i),360);
  if thetap(i)>90
      thetap(i)=180-thetap(i);
  end
  %
  times = sphRes(1,2);
  mults = sphRes(1,3);
  x0(i) = sphRes(1,4);
  y0(i) = sphRes(1,5);
  z0(i) = sphRes(1,6);
  chi2s(i) = sphRes(1,8);
  if z0(i)<REFALT
      z0(i)=REFALT+(REFALT-z0(i));
  end;
  [rhos(i) thetas(i) phis(i)]=Convert2Sph(x0(i),y0(i),z0(i));
end

txtTable = [rhoi*ones(ntries,1) rhos' thetai' thetap' thetas' phii' phip' phis' chi2p' chi2s']

%% Write to file
filename = 'dummyRecons.txt';
fid = fopen( filename, 'w' );        
for l = 1:size( txtTable, 1 )
    fprintf( fid, '%20f', txtTable( l, 1 ) );   % Rho ini
    fprintf( fid, '%20f ',  txtTable( l, 2) );    % Theta ini 
    fprintf( fid, '%8.3f ',  txtTable( l, 3) );    % Theta plan
    fprintf( fid, '%8.3f',  txtTable( l, 4) );    % Theta sph
    fprintf( fid, '%8.3f ', txtTable( l, 5 ) );  % Phi ini
    fprintf( fid, '%8.3f ', txtTable( l, 6 ) );  % Phi plan
    fprintf( fid, '%8.3f ', txtTable( l, 7 ) );  % Phi sph
    fprintf( fid, '%8.3f ', txtTable( l, 8 ) );  %  Chi2 plan
    fprintf( fid, '%8.3f ', txtTable( l, 9 ) );  % Chi2 sph
    fprintf( fid, '%8.3f', txtTable( l, 10 ) );   % 
    fprintf( fid, '\n' );
end
fclose(fid)
