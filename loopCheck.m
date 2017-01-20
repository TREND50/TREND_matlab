function loopCheck()
% OMH 09/03/2015

SharedGlobals;
totEv = zeros(2,length(ALLDETS));

for i = 3562:3733
disp(sprintf('Checking run %d...',i))
  totEv = totEv + checkNbEvent(i);
  totEv'
end
%
totEv'
%
figure(1)
subplot(2,1,1)
semilogy(ALLDETS,totEv(1,:),'b+','MarkerSize',6)
hold on
grid on
semilogy(ALLDETS,totEv(1,:),'k+','MarkerSize',6)
xlabel('AntennaId', labelOpts{:})
ylabel('Nb of events', labelOpts{:})
%
subplot(2,1,2)
d = totEv(1,:)-totEv(2,:);
plot(ALLDETS,d,'r+','MarkerSize',6)
grid on
xlabel('AntennaId', labelOpts{:})
ylabel('\Delta Nb of events', labelOpts{:})
