function [] = eff
% Fit TREND DataChallenge efficiency curve
% SLC+OM
% 14/05/2017

SharedGlobals;

ndays = 80
err_coef = 1

err = 0;
%% Exposure data
x=[0:0.01:30]*1e17;  % E in eV
xs = x*err_coef;
% True case
E=[0 5e16 7e16 8e16 1e17 2e17 3e17 5e17 7e17 1e18 2e18 3e18];
expo=[0 0 5823 17576 38508 262981 611041 1180528 1925580 1956042 2432670 2227544];
expoerr=[0 0 3360 8000 6000 90000 160000 260000 350000 300000 350000 350000];
% E=[0 5e16 8e16 1e17 3e17 5e17 1e18 2e18 3e18];   % 2e18 added by hand for better fit
% expo=[0 0 17576 38508 611041 1180528 1956042 2090000 2227544];
% expoerr=[0 0 8000 6000 160000 100000 300000 350000 400000];
if err == -1
    expo = expo-expoerr;
elseif err == +1
    expo = expo+expoerr;
end
y = interp1( E, expo, xs, 'linear' );
% Ideal cuts
Ecuts=[0 5e16 7e16 8e16 1e17 2e17 3e17 5e17 7e17 1e18 3e18];
expocuts=[0 0 0 35297 72454 551345 1235566 2628733 3722948 4253440 6325718];
expocutserr=[0 0 0 14000 22000 120000 220000 400000 520000 450000 630000];
if err == -1
    expocuts = expocuts-expocutserr;
elseif err == +1
    expocuts = expocuts+expocutserr;
end
%Ecuts=[0 5e16 8e16 1e17 3e17 5e17 1e18 3e18];   % 2e18 added by hand for better fit
%expocuts=[0 0 35297 72454 1235566 2628733 4253440  6325718];
%expocutserr=[0 0 8000 6000 160000 100000 300000 400000];  % Copied from true case
ycuts = interp1( Ecuts, expocuts, xs, 'linear' );
% Ideal detector
Eideal=[0 5e16 7e16 8e16 1e17 2e17 3e17 5e17 7e17 1e18 3e18];
expoideal=[0 0 64930 388565 829610 2675274 5210356 7457073 8758385 9577639 13462189];
expoidealerr=[0 0 16000 50000 80000 300000 400000 700000 900000 700000 1100000];
if err == -1
    expoideal = expoideal-expoidealerr;
elseif err == +1
    expoideal = expoideal+expoidealerr;
end
% Eideal=[0 5e16 8e16 1e17 2e17 3e17 5e17 7e17 1e18 2e18 3e18]; % 2e18 added by hand for better fit
% expoideal=[0 64930 388565 829610 2675274 5210356 7457073 8758385 9577639 11500000 13462189 ];
% expoidealerr=[0 16000 50000 80000 300000 400000 700000 900000 700000 900000 1100000];
yideal = interp1( Eideal, expoideal, xs, 'linear' );

% Erf fit
[ErfPars Res] = FitErf(E/1e16,expo,[1e6,4.3,2.8],0);
ErfPars(1)
ErfPars(2)
ErfPars(3)
yerf = ErfPars(1)*(1+erf((x*1e-16-ErfPars(2))/ErfPars(3)));
yerf = 0.95e6*(1+erf((x*1e-17-4.3)/2.8));
% By hand
yerf = 0.95e6*(1+erf((xs*1e-17-4.3)/2.2));
yerf = 0.95e6*(1+erf((xs*1e-17-4.3)/2.3));
% Sandra's fit
y=0.98e6*(1+erf((x-3.58)/1.66));
y = yerf;

rcuts = expocuts./expoideal;
%r = expo./expoideal;

%% Event number computation
sel_LE = x>=7e16 & x<1e17;
sel_HE = x>=1e17 & x<=3e18;
xg = x/1e9;   % GeV
fluxbef=3.6e6*(xg(sel_LE)).^-3.;   %E<10^17eV
flux=2.3*10^8.6*(xg(sel_HE)).^-3.3;    %E>10^17
nevts(sel_LE) = y(sel_LE).*fluxbef*24*3600*ndays;
nevts(sel_HE) = y(sel_HE).*flux*24*3600*ndays;
nevtscuts(sel_LE) = ycuts(sel_LE).*fluxbef*24*3600*ndays;
nevtscuts(sel_HE) = ycuts(sel_HE).*flux*24*3600*ndays;
nevtsideal(sel_LE) = yideal(sel_LE).*fluxbef*24*3600*ndays;
nevtsideal(sel_HE) = yideal(sel_HE).*flux*24*3600*ndays;
%
totstat = trapz(xg,nevts);
disp(sprintf('Real case: %3.1f events expected in TREND in %d days', totstat,ndays))
totstat = trapz(xg,nevtscuts);
disp(sprintf('Ideal cuts: %3.1f events expected in TREND in %d days', totstat,ndays))
totstat = trapz(xg,nevtsideal);
disp(sprintf('Ideal detector: %3.1f events expected in TREND in %d days', totstat,ndays))

ecuts = [7 9 15 25 40 60 85 200 300]*1e16;
w=zeros(1,length(ecuts)-1);
for i=1:length(ecuts)-1
    sel = x>=ecuts(i) & x<ecuts(i+1);
    w(i)=trapz(xg(sel),nevts(sel));
end
trapz(xg(x>=7e16),nevts(x>=7e16))
w
sum(w)

% weight1e17=55.2543 %between 9e16 and 1.5e17
% weight2e17=67.9120 %between 1.5e17 and 2.5e17
% weight3e17=54.8155 %between 2.5e17 and 4e17
% weight5e17=30.9244 %between 4e17 and 6e17
% weight7e17=18.4067 %between 6e17 and 8.5e18
% weight1e18=13.2683 %between 8.5e17 and 2e18
% weight3e18=1.5305 %between 2e18 and 4e18

%% Plots
figure(1)
subplot(1,1,1)
errorbar(E/1e18/err_coef,expo,expoerr,'+-b','LineWidth',2)
%plot(E/1e18/err_coef,expo,'+-b','LineWidth',3)
grid on
hold on
xlim([0.05 3.1])
errorbar(Eideal/1e18/err_coef,expoideal,expoidealerr,'+-g','LineWidth',2)
errorbar(Ecuts/1e18/err_coef,expocuts,expocutserr,'+-k','LineWidth',2)
%plot(Eideal/1e18/err_coef,expoideal,'+-g','LineWidth',3)
%plot(Ecuts/1e18/err_coef,expocuts,'+-k','LineWidth',3)
plot(x/1e18,ycuts,'r--','LineWidth',1)
plot(x/1e18,yideal,'r--','LineWidth',1)
plot(x/1e18,y,'r--','LineWidth',1)
plot(x/1e18,yerf,'m','LineWidth',3)

xlabel('Energy [EeV]',labelOpts{:})
ylabel('Apperture [m^2.sr]',labelOpts{:})
%title('TREND exposure', labelOpts{:})
% subplot(2,1,2)
% semilogx(E/1e18,expo,'+-b','LineWidth',2)
% grid on
% hold on
% xlim([0.05 3.1])
% semilogx(Ecuts/1e18,expocuts,'+-k','LineWidth',2)
% semilogx(Eideal/1e18,expoideal,'+-g','LineWidth',2)
% %semilogx(x/1e18,y,'r','LineWidth',2)
% xlabel('Energy [EeV]',labelOpts{:})
% ylabel('Apperture [m^2.sr]',labelOpts{:})

% figure(2)
% semilogx(E/1e18/err_coef,r,'+-b','LineWidth',2)
% grid on
% hold on
% semilogx(E/1e18/err_coef,rcuts,'+-k','LineWidth',2)
% xlim([0.05 3.1])
% xlabel('Energy [EeV]',labelOpts{:})
% ylabel('Apperture ratio to ideal',labelOpts{:})

% figure(3)
% subplot(2,1,1)
% errorbar(Eideal/1e18,expoideal,expoidealerr,'+-k','LineWidth',2)
% grid on
% hold on
% xlim([0.05 3.1])
% plot(x/1e18,y,'r','LineWidth',2)
% xlabel('Energy [EeV]',labelOpts{:})
% ylabel('Apperture [m^2.sr]',labelOpts{:})
% %title('Ideal case', labelOpts{:})
% subplot(2,1,2)
% semilogx(Eideal/1e18,expoideal,'+-k','LineWidth',2)
% grid on
% hold on
% xlim([0.05 3.1])
% semilogx(x/1e18,y,'r','LineWidth',2)
% xlabel('Energy [EeV]',labelOpts{:})
% ylabel('Apperture [m^2.sr]',labelOpts{:})

figure(7)
loglog(xg(sel_LE),fluxbef,'k','LineWidth',2)
hold on
grid on
loglog(xg(sel_HE),flux,'k','LineWidth',2)
xlim([5e7 3e9])
xlabel('Energy [GeV]',labelOpts{:})
ylabel('dN/dEdtd\OmegadS [GeV^{-1}s^{-1}m^{-2}sr^{-1}]',labelOpts{:})
title('CR flux', labelOpts{:})

figure(11)
loglog(x(sel_LE),nevts(sel_LE),'b','LineWidth',2)
hold on
grid on
loglog(x(sel_HE),nevts(sel_HE),'b','LineWidth',2)
loglog(x(sel_LE),nevtscuts(sel_LE),'k','LineWidth',2)
loglog(x(sel_HE),nevtscuts(sel_HE),'k','LineWidth',2)
loglog(x(sel_LE),nevtsideal(sel_LE),'g','LineWidth',2)
loglog(x(sel_HE),nevtsideal(sel_HE),'g','LineWidth',2)
xlabel('Energy [eV]',labelOpts{:})
ylabel('dN/dE [GeV^{-1}]',labelOpts{:})
title(sprintf('Nb of events in TREND in %d days',ndays), labelOpts{:})
xlim([5e16 3.1e18])
