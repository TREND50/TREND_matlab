function[]=cutseffect
% Detail of effect of different cuts on EAS selection
% SLC 29/05/2017

SharedGlobals;

%1 mutiplicité>4
%2 rayon>500
%3 chi2<30
%4 thetap<81
%5 mean(barydist)>500
%6 amplitude
%7 status, bads<=1
%8 ncomdir <=3
%9 ncommon <=10

%10 bads<=0
%11 chi2 <=30
%12 radius >=3000
%13 ratio amplitude >=1
%14 theta<=81 et theta>0
%15 dirneib
%16 neib
%(17 if period>=9 ...)

col = {'k','b','m','r','g','y','k--','b--','m--','r--','g--','y--','k.-','b.-','m.-','r.-','g.-','y.-',}
En=[8e16 1e17 2e17 3e17 5e17 7e17 1e18 3e18];

% Results of simulations
a=[9973 35 18 31 27 27 12 12 12 12 12 12 12 12 12 12 12 10 10 10 9 9 9;   
9746 65 37 58 51 50 26 25 25 25 25 25 24 24 24 24 24 22 22 22 22 20 20;   
3989 113 81 106 88 86 58 58 57 57 57 57 55 55 55 54 54 52 52 51 43 38 38;  
3341 204 148 192 155 150 107 105 103 103 103 103 99 99 99 98 98 95 95 95 83 74 74;
2796 228 187 206 177 170 138 137 136 132 128 128 120 120 120 118 118 113 113 113 94 87 87;
1981 234 203 230 201 198 164 163 163 161 156 156 143 141 140 135 135 132 132 132 114 105 105;
3250 439 372 425 366 331 286 281 276 271 255 255 232 232 232 228 228 222 222 221 184 169 169;
2899 544 473 536 461 425 363 361 339 329 302 282 263 263 262 213 213 203 203 203 176 166 166]

ref=a(:,2);  
for i=1:length(En)
  r(i,:)=a(i,:)/ref(i);
end

figure(37)
semilogx(En,r(:,4:end),'-*','Linewidth',2)
xlim([7e16 4e18])
grid on

xlabel('Energy [eV]',labelOpts{:})
ylabel('N_{sel}/NOnlineCoincs_{4ants}',labelOpts{:})

%cuts 1 a 17:
r(:,7:end)


