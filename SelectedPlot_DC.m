function [] = SelectedPlot_DC(E)

SharedGlobals
SIMCODE='ZHAIRESall'

selfile = load([CAND_PATH E '/log_candanalysisend.txt'], 'r' );
theta=selfile(:,2);
phi=mod(selfile(:,3)-90,360);  % use TREND convention
sel=selfile(:,33);
thetarecons=selfile(:,15);
phirecons=selfile(:,16);
ntrigs=selfile(:,8);
nindst=selfile(:,14);
baryOK=selfile(:,21); % Surviving Bary cut
statusOK=selfile(:,22);  %Status
patternOK=selfile(:,23);  %Pattern
Chi2OK=selfile(:,27);  % Chi2plan
RcutOK = selfile(:,28);  %Rcut
ThetaOK=selfile(:,30);  % Theta
DirNeiOK=selfile(:,31); %DirNeighbourgh
NeiOK=selfile(:,32);  % Neighbourgh

recoverfile=load([DST_PATH E '/log_recoverall.txt'], 'r' );
thetabef=recoverfile(:,2);
phibef=recoverfile(:,3);
nantscoinc=recoverfile(:,10);
nantspp=recoverfile(:,12);
trigbef=recoverfile(:,8);

nsimu=length(thetabef);
n4ants=length(trigbef(trigbef>=4));
n5ants=length(trigbef(trigbef>=5));
ncoinc=length(nantscoinc(nantscoinc~=0))+length(nindst)
npp=length(nantspp(nantspp~=0))+length(nindst)
ndst=length(nindst)
nBaryOK=sum(baryOK)
nStatusOK=sum(statusOK)
nPatternOK=sum(patternOK)
nRcutOK=sum(RcutOK)
nChi2OK=sum(Chi2OK)
nThetaOK=sum(ThetaOK)
nDirNeiOK=sum(DirNeiOK)
nNeiOK=sum(NeiOK)
surv=sum(sel)

display('nsimu,n4ants,n5ants,ncoinc,npp,ndst:')
display(sprintf('%i %i %i %i %i %i',nsimu,n4ants,n5ants,ncoinc,npp,ndst))
display('Surviving Bary,Status,Pattern,Chi2,Rcut,Theta,DirNei,Nei,Sel:')
display(sprintf('%i %i %i %i %i %i %i %i %i',nBaryOK,nStatusOK,nPatternOK,nChi2OK,nRcutOK,nThetaOK,nDirNeiOK,nNeiOK,surv))

figure(2)
PrepareSkyPlot(2);
polar( phibef*pi/180, thetabef, 'ko' );
fig=polar( phi(sel==1)*pi/180, theta(sel==1), 'go');
set(fig, 'MarkerFaceColor','g' );  
saveas(fig,sprintf('cible%s%s.png',SIMCODE,E))
hold off

figure(1)
hist(phi(sel==1),9)
saveas(figure(1),sprintf('cut9phi%s%s.png',SIMCODE,E))

figure(3)
hist(phi,9)
saveas(figure(3),sprintf('coincphi%s%s.png',SIMCODE,E))

figure(7)
hist(theta(sel==1),9)
saveas(figure(7),sprintf('cut9theta%s%s.png',SIMCODE,E))

figure(8)
hist(theta,9)
saveas(figure(8),sprintf('coinctheta%s%s.png',SIMCODE,E))

figure(4)
plot(theta(sel==1),thetarecons(sel==1),'g*')
hold on
plot(theta(sel==0),thetarecons(sel==0),'k*')
x=[0:1:90];
plot(x,x,'k')
grid on
xlabel('Theta_{True}', labelOpts{:})
ylabel('Theta_{Recons}', labelOpts{:})
saveas(figure(4),sprintf('thetathetarecons%s%s.png',SIMCODE,E))

figure(5)
hist(nindst./ntrigs,10);
saveas(figure(5),sprintf('nantsratio%s%s.png',SIMCODE,E))

figure(6)
plot(phi(sel==1),mod(phirecons(sel==1),360),'g*')
hold on
plot(phi(sel==0),mod(phirecons(sel==0),360),'k*')
x=[0:1:360];
plot(x,x,'k')
grid on
xlabel('Phi_{True} [deg]', labelOpts{:})
ylabel('Phi_{recons} [deg]', labelOpts{:})
saveas(figure(6),sprintf('phiphirecons%s%s.png',SIMCODE,E))

