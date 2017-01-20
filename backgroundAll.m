function []=backgroundAll()
% Background sensitivity plots
% OMH 14/06/2011

SharedGlobals;
antennas = [101:138 140 148:158];
nruns = 3050:3053;
%nruns = 2693:3086;

cod = {'+k','+b','+m','+r','+g','^k','^b','^m','^r','^g','vk','vb','vm','vr','vg','pk','pb','pm','pr','pg','hk','hb','hm','hr','hg',
       '+k','+b','+m','+r','+g','^k','^b','^m','^r','^g','vk','vb','vm','vr','vg','pk','pb','pm','pr','pg','hk','hb','hm','hr','hg'};

for i = 1:length(antennas)
   disp(sprintf('backgroundAll: performing full analysis for antenna %d',antennas(i)))
   [sig,t,lst]=stdbackground(nruns,antennas(i));
   if size(sig,1)==0
       continue
   end
   figure(i)
   set(i,'Name',sprintf('Background - Antenna %d',antennas(i)),'NumberTitle','off','Position',[1 1 scrsz(3)/2 scrsz(4)]);
   subplot(2,1,1)
   plot(lst,sig,'+k')
   xlabel('LST [h]',labelOpts{:})
   ylabel('Sigma [LSB]',labelOpts{:})
   axis([0, 24, 0, 12])
   grid on
   subplot(2,1,2)
   plot(t,sig,'+k')
   datetick
   ylim([0 12])
   xlabel('Time [h]',labelOpts{:})
   ylabel('Sigma [LSB]',labelOpts{:})
   grid on
   %pause
   close(i)
   %
   figure(1)
   set(1,'Name','Background all','NumberTitle','off')
   plot(lst,sig,cod{i})
   hold on
   grid on
   xlabel('LST [h]',labelOpts{:})
   ylabel('Sigma [LSB]',labelOpts{:})
end