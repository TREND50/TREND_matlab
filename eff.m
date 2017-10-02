function [] = eff
% Fit TREND DataChallenge efficiency curve
% SLC+OM
% 14/05/2017

SharedGlobals;

ndays = 80
err_coef = 0.8

err = 0;
%% Exposure data
x=[0:0.01:30]*1e17;  % E in eV
xs = x*err_coef;
% Proton Period 6 - True case
a(1,:) = [7e16, 5754.9777848744006, 2432.9518469638338, 9077.0037227849662];
a(end+1,:) = [8e16, 17173.830251649244, 9415.8505082187276, 24931.809995079762];
a(end+1,:) = [1e17, 38987.00461787219, 23996.961505230654, 53977.047730513732];
a(end+1,:) = [2e17, 260886.89065799001, 174754.92830298044, 347018.85301299958];
a(end+1,:) = [3e17, 611191.25119069754, 465755.29919609026, 756627.20318530477];
a(end+1,:) = [5e17, 1156636.0479507558, 924518.59068996971, 1388753.5052115421];
a(end+1,:) = [7e17, 2038158.072360653, 1651914.8697151789, 2424401.2750061275];
a(end+1,:) = [1e18, 1960580.3642857312, 1667762.6122867302, 2253398.1162847327];
a(end+1,:) = [2e18, 2361387.2826766241, 2022140.764262483, 2700633.8010907657];
a(end+1,:) = [3e18, 2171572.1715591662, 1855949.8067556992, 2487194.5363626331];
E = a(:,1);
expo = a(:,2);
expoerrsup = a(:,4)-a(:,2);
expoerrinf = a(:,2)-a(:,3);
expoerr = (expoerrsup+expoerrinf)/2;
y = interp1( E, expo, xs, 'linear' );


% Fe period 6
a = [];
a(1,:) = [1e17,25323.960185008025, 13328.312997147412, 37319.607372868639];
a(end+1,:) = [3e17,669371.37768418854, 534545.22682648245, 804197.52854189475];
a(end+1,:) = [5e17,1228751.0237451992, 975525.29059844417, 1481976.7568919547];
a(end+1,:) = [1e18,2517456.3801536262, 2155778.7296749926, 2879134.0306322593];
EFe = a(:,1);
expoFe = a(:,2);
expoFeerrsup = a(:,4)-a(:,2);
expoFeerrinf = a(:,2)-a(:,3);
expoFeerr = (expoFeerrsup+expoFeerrinf)/2;

% Ideal cuts Period 6
a = [];
a(1,:) = [7e16, 7714.0335069772636, 3857.9691033188255, 11570.097910635703];
a(end+1,:) = [8e16, 35050.831362124125, 21083.144179519637, 49018.51854472862];
a(end+1,:) = [1e17,71840.834881414834, 49897.971582129794, 93783.698180699866];
a(end+1,:) = [2e17,543342.88451279502, 423166.42953207035, 663519.33949351986];
a(end+1,:) = [3e17,1226459.6438238451, 1012024.2665327355, 1440895.0211149547];
a(end+1,:) = [5e17,2618751.023764574, 2225688.9377190717, 3011813.1098100767];
a(end+1,:) = [7e17,3846809.5924892053, 3326930.3469126327, 4366688.8380657788];
a(end+1,:) = [1e18,4194901.8854535026, 3761073.9206238505, 4628729.8502831543];
a(end+1,:) = [2e18, 5780553.9111783858, 5208891.8201225121, 6352216.0022342615];
a(end+1,:) = [3e18, 6442186.6864160709, 5820142.6518890038, 7064230.7209431389];
Enocuts = a(:,1);
exponocuts = a(:,2);
exponocutserrsup = a(:,4)-a(:,2);
exponocutserrinf = a(:,2)-a(:,3);
exponocutserr = (exponocutserrsup+exponocutserrinf)/2;
if err == -1
    expocuts = expocuts-expocutserr;
elseif err == +1
    expocuts = expocuts+expocutserr;
end
ynocuts = interp1( Enocuts, exponocuts, xs, 'linear' );
ynocuts(isnan(ynocuts)) = 0;

% Ideal detector Period 6
%a(1,:) = [5e16,65186.579439257825, 48519.081861045313, 81854.077017470307];
a = [];
a(1,:) = [7e16, 275812.89135843155, 238010.01762385212, 313615.76509301097];
a(end+1,:) = [8e16, 389023.61257966381, 338389.67349632649, 439657.55166300107];
a(end+1,:) = [1e17, 821111.96290013869, 742338.87393297395, 899885.05186730332];
a(end+1,:) = [2e17, 2611884.4746014378, 2333298.35665621, 2890470.5925466651];
a(end+1,:) = [3e17,5311879.0752844503, 4875295.9128011372, 5748462.2377677634];
a(end+1,:) = [5e17,7517272.4278993644, 6859840.2682627309, 8174704.5875359979];
a(end+1,:) = [7e17,8953225.7888710164, 8065501.7045939304, 9840949.8731481023];
a(end+1,:) = [1e18,9741657.2999682277, 9042783.2689351328, 10440531.331001326];
a(end+1,:) = [2e18,12036320.781420546, 11188756.633816075, 12883884.929025019];
a(end+1,:) = [3e18,14256592.443819441, 13219180.272674572, 15294004.614964306];
Eideal = a(:,1);
expoideal = a(:,2);
expoidealerrsup = a(:,4)-a(:,2);
expoidealerrinf = a(:,2)-a(:,3);
expoidealerr = (expoidealerrsup+expoidealerrinf)/2;
e = load('../soft_cal/e_ideal.txt');
f = load('../soft_cal/fit_ideal.txt');
yideal = interp1( e, f, xs, 'linear' );
if err == -1
    expoideal = expoideal-expoidealerr;
elseif err == +1
    expoideal = expoideal+expoidealerr;
end
yideal(isnan(yideal)) = 0;

% p Period 8
a = [];
a(1,:) = [7e16, 0,0,0];                                                              % Fake
a(end+1,:) = [8e16, 1946.2310267221776, 0.12094776923228573, 3892.3411056751229];
a(end+1,:) = [1e17,15884.838356062608, 7145.8961858448529, 24623.780526280359];
a(end+1,:) = [2e17,171502.73049910343, 112751.42670327882, 230254.03429492805];       
a(end+1,:) = [3e17,438081.23851231206, 326843.84826787992, 549318.62875674432];
a(end+1,:) = [5e17,1055501.3356437786, 838606.98575439188, 1272395.6855331652];
a(end+1,:) = [7e17,1223357.9197879077, 896197.9944760364, 1550517.8450997788];
a(end+1,:) = [1e18,1220241.1590411672, 989297.57922642096, 1451184.7388559133];
a(end+1,:) = [3e18,1220241.1590411672, 989297.57922642096, 1451184.7388559133];      % Fake

E8 = a(:,1);
expo8 = a(:,2);
expoerrsup8 = a(:,4)-a(:,2);
expoerrinf8 = a(:,2)-a(:,3);
expo8err = (expoerrsup8+expoerrinf8)/2;
y8 = interp1( E8, expo8, xs, 'linear' );
y8(isnan(y8)) = 0;

%% Erf fit
% [ErfPars Res] = FitErf(E/1e16,expo,[1e6,4.3,2.8],0);
% yerf = ErfPars(1)*(1+erf((x*1e-16-ErfPars(2))/ErfPars(3)));
% yerf = 0.95e6*(1+erf((x*1e-17-4.3)/2.8));
% Sandra's fit
A = 0.9800;
B = 3.5800;
C = 1.6600;
yerf = A*1e6*(1+erf((xs*1e-17-B)/C));

%% Valentin's fit
e = load('../soft_cal/e_norm.txt');
f = load('../soft_cal/fit_norm.txt');
yf = interp1( e, f, xs, 'linear' );


%y = yerf;
y = yf;
y(isnan(y)) = 0;

rcuts = exponocuts./expoideal;
r = expo./expoideal;

%% Event number computation
sel_LE = x>=7e16 & x<1e17;
sel_HE = x>=1e17 & x<=3e18;
xg = x/1e9;   % GeV
fluxbef=3.6e6*(xg(sel_LE)).^-3.;   %E<10^17eV
flux=2.3*10^8.6*(xg(sel_HE)).^-3.3;    %E>10^17
nevts(sel_LE) = y(sel_LE).*fluxbef*24*3600*ndays;
nevts(sel_HE) = y(sel_HE).*flux*24*3600*ndays;
nevtsnocuts(sel_LE) = ynocuts(sel_LE).*fluxbef*24*3600*ndays;
nevtsnocuts(sel_HE) = ynocuts(sel_HE).*flux*24*3600*ndays;
nevtsideal(sel_LE) = yideal(sel_LE).*fluxbef*24*3600*ndays;
nevtsideal(sel_HE) = yideal(sel_HE).*flux*24*3600*ndays;
%
ndaysP8 = 28.6;
nevts8(sel_HE) = y8(sel_HE).*flux*24*3600*ndaysP8;
nevts8(sel_LE) = y8(sel_LE).*fluxbef*24*3600*ndaysP8;
totstatP8 = trapz(xg,nevts8);
%
totstat = trapz(xg,nevts);
disp(sprintf('Real case: %3.1f events expected in TREND in %d days', totstat,ndays))
totstats = trapz(xg,nevtsnocuts);
disp(sprintf('Ideal cuts: %3.1f events expected in TREND in %d days', totstats,ndays))
totstatd = trapz(xg,nevtsideal);
disp(sprintf('Ideal detector: %3.1f events expected in TREND in %d days', totstatd,ndays))
%
%totstat./totstats
%totstats./totstatd

totstat = trapz(xg,nevts);
disp(sprintf('Real case: %3.1f events expected in TREND in %3.1f days of period 8', totstatP8,ndaysP8))

%ecuts = [7 9 15 25 40 60 85 200 300]*1e16;
d = [1; diff(E)]/2;
ecuts = [E-d; E(end)+100];
w=zeros(1,length(ecuts)-1);
for i=1:length(ecuts)-1
    sel = x>=ecuts(i) & x<ecuts(i+1);
    w(i)=trapz(xg(sel),nevts(sel));
end
trapz(xg(x>=7e16),nevts(x>=7e16))
w
sum(w)


%% Plots
figure(1)
subplot(2,1,1)
errorbar(E/1e18/err_coef,expo,expoerr,'+-b','LineWidth',2)
%plot(E/1e18/err_coef,expo,'+-b','LineWidth',3)
grid on
hold on
xlim([0.05 3.1])
errorbar(Eideal/1e18/err_coef,expoideal,expoidealerr,'+-g','LineWidth',2)
errorbar(Enocuts/1e18/err_coef,exponocuts,exponocutserr,'+-m','LineWidth',2)
%errorbar(EFe/1e18/err_coef,expoFe,expoFeerr,'+k','LineWidth',2)
%errorbar(E8/1e18/err_coef,expo8,expo8err,'+r','LineWidth',2)
%plot(Eideal/1e18/err_coef,expoideal,'+-g','LineWidth',3)
%plot(Ecuts/1e18/err_coef,expocuts,'+-k','LineWidth',3)
%plot(x/1e18,ynocuts,'r--','LineWidth',1)
%plot(x/1e18,yideal,'r--','LineWidth',1)
%plot(x/1e18,y,'r--','LineWidth',1)
%plot(x/1e18,yerf,'b','LineWidth',3)
%plot(x/1e18,y8,'r','LineWidth',3)
plot(xs/1e18,yf)
xlabel('Energy [EeV]',labelOpts{:})
ylabel('Apperture [m^2.sr]',labelOpts{:})
%title('TREND apperture', labelOpts{:})
subplot(2,1,2)
loglog(E/1e18/err_coef,expo,'*b','LineWidth',2,'MarkerSize',12)
grid on
hold on
xlim([0.05 3.1])
loglog(Eideal/1e18,expoideal,'*-g','LineWidth',2,'MarkerSize',12)
loglog(Enocuts/1e18,exponocuts,'*-m','LineWidth',2,'MarkerSize',12)
loglog(EFe/1e18/err_coef,expoFe,'*k','LineWidth',2,'MarkerSize',12)
loglog(E8/1e18/err_coef,expo8,'*r','LineWidth',2,'MarkerSize',12)
%loglog(x/1e18,yerf,'b','LineWidth',3,'MarkerSize',8)
loglog(xs/1e18/err_coef,yf,'b','LineWidth',3,'MarkerSize',8)
%semilogx(x/1e18,y,'r','LineWidth',2)
xlabel('Energy [EeV]',labelOpts{:})
ylabel('Apperture [m^2.sr]',labelOpts{:})

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
%loglog(x(sel_LE),nevtsnocuts(sel_LE),'k','LineWidth',2)
%loglog(x(sel_HE),nevtsnocuts(sel_HE),'k','LineWidth',2)
%loglog(x(sel_LE),nevtsideal(sel_LE),'g','LineWidth',2)
%loglog(x(sel_HE),nevtsideal(sel_HE),'g','LineWidth',2)
xlabel('Energy [eV]',labelOpts{:})
ylabel('dN/dE [GeV^{-1}]',labelOpts{:})
title(sprintf('Nb of events in TREND in %d days',ndays), labelOpts{:})
xlim([5e16 3.1e18])
