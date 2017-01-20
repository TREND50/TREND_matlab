function [pos_ant,podN,pos_pod,cable,delay,sciID,logID]=SetupCharacteristics2010(Antennas)
% Caractéristiques du réseau d'antennes en fonction du run étudié
% Réanalyse des runs hybrides 2010 (R1788-2097)
% OMH 01/09/2011

disp '************ Warning!!! Using SetupCharacteristics2010! ***************** '
pause(5)

for i=1:length(Antennas)
    % Position: coord_crosspoint3.txt
    % Delays antennas: delays_crosspoint.txt
    % Delays scints: delays_scintcor.txt (nrun<1900)
    %                delays_scint2cor.txt (nrun>1900)
    switch Antennas(i)
        %% Antennas            
        case 136
            pos_ant(i,:)=[333.0 -14.1 2638.6];
            podN(i) = 136;
            pos_pod(i,:)=[349.0       0.0	 2636.2];
            cable(i)=0;
            delay(i)=-187.4;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 137
            pos_ant(i,:)=[241.0 30.7 2637.6 ];
            podN(i) = 137;
            pos_pod(i,:)=[269.0       0.0	 2635.0];
            cable(i)=0;
            delay(i)=-193.9;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 138
            pos_ant(i,:)=[233.8 -40.7 2638.4 ];
            podN(i) = 138;
            pos_pod(i,:)=[249.0       0.0	 2634.7];
            cable(i)=0;
            delay(i)=-164.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 140
            pos_ant(i,:)=[109.0 -100.0 2639.0];
            podN(i) = 140;
            pos_pod(i,:)=[109.0       0.0     	 2633.1];
            cable(i)=0;
            %delay(i)=0;  
            delay(i)=5.2*0.7;   %run>1992
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 148
            pos_ant(i,:)=[61.4 571.7 2645.1];
            podN(i) = 148;
            pos_pod(i,:)=[0.0	   532.0     	 2638.7];
            cable(i)=0;
            delay(i)=627.8;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 149
            pos_ant(i,:)=[-18.1 537.1 2644.6];
            podN(i) = 149;
            pos_pod(i,:)=[0.0	   502.7     	 2639.0];
            cable(i)=0;
            delay(i)=576.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 150
            pos_ant(i,:)=[100.3 472.7 2641.9];
            podN(i) = 150;
            pos_pod(i,:)=[0.0	   473.3     	 2639.3];
            cable(i)=0;
            delay(i)=563.4+7.6*0.7; %nrun>1700
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 151
            pos_ant(i,:)=[35.0 462.0 2642.0];
            podN(i) = 151;
            pos_pod(i,:)=[0.0	   443.9     	 2639.6 ];
            cable(i)=0;
            delay(i)=490.8;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 152
            pos_ant(i,:)=[-22.7 399.5 2641.9];
            podN(i) = 152;
            pos_pod(i,:)=[0.0	   414.6     	 2639.9];
            cable(i)=0;
            delay(i)=462.8;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 153
            pos_ant(i,:)=[92.0 354.5 2638.9];
            podN(i) = 153;
            pos_pod(i,:)=[0.0	   385.2     	 2639.3];
            cable(i)=0;
            delay(i)=474.8;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 154
            pos_ant(i,:)=[-24.8 268.9 2634.9];
            podN(i) = 154;
            pos_pod(i,:)=[0.0	   238.4     	 2632.7];
            cable(i)=0;
            delay(i)=287.9;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 155
            pos_ant(i,:)=[66.7 201.6 2635.3];
            podN(i) = 155;
            pos_pod(i,:)=[0.0	   209.0     	 2631.4];
            cable(i)=0;
            delay(i)=269.6+5.4*0.7; %nrun>1700
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 156
            pos_ant(i,:)=[-11.3 11.3 2631.5];
            podN(i) = 156;
            pos_pod(i,:)=[0.0	  -055.3     	 2633.4];
            cable(i)=0;
            %delay(i)=134.6;
            delay(i)=134.6+6.0*0.7; %nrun>1992
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 157
            pos_ant(i,:)=[31.4 -89.1 2634.9];
            podN(i) = 157;
            pos_pod(i,:)=[0.0	  -172.8     	 2636.3];
            cable(i)=0;
            %delay(i)=250.1;
            delay(i)=250.1+6.0*0.7; %nrun>1992;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 158
            pos_ant(i,:)=[81.1 -186.7 2637.7];
            podN(i) = 158;
            pos_pod(i,:)=[0.0	  -202.1     	 2637.1];
            cable(i)=0;
            delay(i)=298.0;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1;     
        
        %% Scintillators... recompute delays
        case 139  
            pos_ant(i,:)=[149.0  -30.0 2633.5];
            podN(i) = 139;
            pos_pod(i,:)=[149.0       0.0	 2633.5];
            cable(i)=0;
            delay(i)=-32.2; %run<1900
            delaycorr(i)=0;
            %delay(i)=-71.4;  %run>1900
            sciID(i)=1;
            logID(i)=0; 
            
            
        case 111
            pos_ant(i,:)=[4.5     197.4  2631.4];
            podN(i) = 155;
            pos_pod(i,:)=[0.0	   209.0     	 2631.4];
            cable(i)=0;
            delay(i)=240.4;
            delaycorr(i)=0;
            sciID(i)=1;            
            logID(i)=0; 
            
        case 110
            pos_ant(i,:)=[0.0     -32.2  2633.4];
            podN(i) = 156;
            pos_pod(i,:)=[0.0	  -055.3     	 2633.4];
            cable(i)=0;
            delay(i)=81.0;
            delaycorr(i)=0;
            sciID(i)=1;
            logID(i)=0; 
            
        otherwise
            display(sprintf('SetupCharacteristics error : Unknow antenna %3d',ant(i)))
            break;
    end;
    
end;