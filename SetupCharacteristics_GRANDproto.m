function [pos_det,podN,machines,delay,cable,isSci]=SetupCharacteristics_GRANDproto(dets,nrun)
% Caractéristiques du réseau d'antennes en fonction du run 
% 14/04/10, TS
% Adapted for GRANDproto
% % Last modification : OMH 08/08/2015
% "delay" units = samples = 5ns

SharedGlobals;
%nrun_temp=num2str(nrun);
%nruntrue=str2num(nrun_temp(1:4));

if nrun<6700
  disp(sprintf('SetupCharacteristics_GRANDproto.m not valid for runs below R6700.\nCalling SetupCharacteristics for previous setup.'))
  return
end



for i=1:length(dets)
    
    switch dets(i)
        %% Antennas  
        case 1   % Former B4 ==> now A127 AERA electronics
            pos_det(i,:)=[19.0	441.9	2634.0];
            podN(i,:) = [151, 151, 151];  %S11
            if nrun>7120
                machines(i,:) = [128, 131, 132]; % Since April 2016
            else
                machines(i,:) = [129, 131, 132];  % September 2015
            end
            cable(i,:)= [+14, +14, +12];  % 9.5+9.2 for X, same for Y, 9.5+6.2 for Z
            fiberdelay(i,:) = [-70, -71 , -71]; % 
            isSci(i)= 0;
        
        case 2% C5
            pos_det(i,:) = [152.8,466.3,2633.6];
            %pos_ant(i,:)=[147, 461, 2640];  Prelim meas.
            podN(i,:) = [150, 150, 150];   %S10
            machines(i,:) = [134,135,136];
            cable(i,:) = [122, 124, 121]; %R6899, 6902, 6903
            fiberdelay(i,:) = [-99, -99 , -99]; %R6886, 6888, 6889  
            isSci(i)= 0;
        
        case 3  % Former C3  )==> now A105 AERA 
            pos_det(i,:)=[-66.6	477.7	2636.1];
            podN(i,:) = [149,149,149];
            machines(i,:) = [137,138,140];
            cable(i,:) = [+56,+61,+56]; %R6906, 6907, 6908
            fiberdelay(i,:) = [-122, -124 , -122]; %R6905,6910,6909     
            isSci(i)= 0;
        
        case 4  % Former C4    DEAD
            pos_det(i,:)=[-6.7	512.2	2636.2];
            %pos_ant(i,:)=[-10, 518, 2639];  % Prelim meas.
            podN(i,:) = [148,148,148];
            machines(i,:) = [151,154,155];
            cable(i,:) = [+8,+8,+8]; %Mesure directe: cable = 10m
            fiberdelay(i,:) = [-162, -160, -162]; %R6911-6913                     
            isSci(i)= 0;

        case 5  % Former A5 - 2D antenna - DEAD
            pos_det(i,:)=[-42.0	635.0	2640.1];
            podN(i,:) = [147, 146, 0];
            machines(i,:) = [148, 149, 0];
            cable(i,:) = [44.5, 81, 0];  % R6974, R6967
            fiberdelay(i,:) = [-298, -347.5 , 0]; % R6922, R6938
            isSci(i)= 0;
        
        case 6  % Former A3  ==> Now AERA126
            pos_det(i,:)=[50.9	657.6	2640.3];
            podN(i,:) = [146, 146, 146];
            machines(i,:) = [156, 157, 158];
            cable(i,:) = [72.5, 70, 69.5];  %R6971,R6969,R6973 
            fiberdelay(i,:) = [-347.5, -348.5 , -347.5]; %R6935-7
            isSci(i)= 0;
               
            
        %% Scints
        case 111% S1 (West)
            pos_det(i,:)=[-74.7	391.2	2633.7];
            podN(i,:) = [153,0,0];
            machines(i,:) = [111,0,0];
            cable(i,:) = [+68,0,0]; %R6871
            fiberdelay(i,:) = [+1,0,0]; %R6867
            isSci(i)= 1;
            
        case 112% S2 (Center)
            pos_det(i,:)=[15.4	395.9	2632.9];
            podN(i,:) = [153,0,0];
            machines(i,:) = [112,0,0];
            cable(i,:) = [+32,0,0 ]; %R6872
            fiberdelay(i,:) = [0,0,0]; %R6869   
            isSci(i)= 1;
            
        case 113% S3 (East)  % Green fiber @ pod S13 is reference channel.
            pos_det(i,:)=[105.1	395.7	2632.4];
            podN(i,:) = [153, 0,0];
            machines(i,:) = [113,0,0];
            cable(i,:) = [+80,0,0 ]; %R6874
            fiberdelay(i,:) = [0,0,0]; %Reference   
            isSci(i)= 1;
            
        case 114% S4 (West)
            pos_det(i,:)=[-74.4	585.3	2639.1];
            podN(i,:) = [147,0,0];
            machines(i,:) = [114,0,0];
            cable(i,:) = [+89,0,0]; %R6943
            fiberdelay(i,:) = [-298,0,0]; % R6920
            isSci(i)= 1;
            
        case 115% S5 (Center)
            pos_det(i,:)=[-4.2	595.3	2638.5];
            podN(i,:) = [147,0,0];
            machines(i,:) = [115,0,0];
            cable(i,:) = [+62,0,0 ]; %R6944
            fiberdelay(i,:) = [-297,0,0]; % R6923
            isSci(i)= 1;
            
        case 116% S6 (East)
            pos_det(i,:)=[100.1	581.8	2637.1];
            podN(i,:) = [147, 0,0];
            machines(i,:) = [116,0,0];
            cable(i,:) = [+112,0,0 ];  %R6945
            fiberdelay(i,:) = [-296,0,0]; % R6919   
            isSci(i)= 1;
            
        otherwise
            display(sprintf('SetupCharacteristics error : Unknow detector %3d',dets(i)))
            %break;
    end;  
end;
delay = cable+fiberdelay;
