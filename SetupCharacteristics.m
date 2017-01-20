function [pos_ant,podN,pos_pod,cable,delay,delaycorr,sciID,logID]=SetupCharacteristics(Antennas,nrun)
% Caract�ristiques du r�seau d'antennes en fonction du run �tudi�
% 14/04/10, TS
% Last modification : OMH 08/06/2011
%If the run is sliced, the nrun should 6 numbers,be like 469000,
%19/10/2014, Jianli

SharedGlobals;
nrun_temp=num2str(nrun);
nruntrue=str2num(nrun_temp(1:4));
gps = 1;
for i=1:length(Antennas)
    
    switch Antennas(i)
        %% Antennas
        case 101
            if gps
                pos_ant(i,:)=[1466.7 3.8 2653.0];
            else
                pos_ant(i,:)=[1471.5 1.4 2653.0];
            end
            podN(i) = 101;
            pos_pod(i,:)=[1560.0       0.0     	 2653.7];
            cable(i)=0;  %samples
            delay(i)=0+3*PROPAGT;  % samples
            delaycorr(i)=0;
            % Cable is 6m long on A101 instead of 3m elsewhere.
            sciID(i)=0;  % 0: antenna 1: scint
            logID(i)=0;  % 0: butterfly 1: log periodical
            
        case 102
            if gps
                pos_ant(i,:)=[1464.9 68.2 2653.0];
            else
                pos_ant(i,:)=[1467.0   66.0 2653.0];
            end
            podN(i) = 102;
            pos_pod(i,:)=[1600.0       0.0     	 2654.1];
            cable(i)=0;
            delay(i)=98.7;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 103
            if gps
                pos_ant(i,:)=[1594.0 -205.2 2653.0];
            else
                pos_ant(i,:)=[1599.0 -207.0 2653.0];
            end
            podN(i) = 103;
            pos_pod(i,:)=[1700.0       0.0     	 2655.1];
            cable(i)=0;
            delay(i)=281.7;
            delaycorr(i)=+4.9;
            sciID(i)=0;
            logID(i)=0; 
            
        case 104
            if gps
                pos_ant(i,:)=[1592.5 194.8 2656.0];
            else
                pos_ant(i,:)=[1597.0  192.0 2656.0];
            end
            podN(i) = 104;
            pos_pod(i,:)=[1720.0       0.0     	 2655.3];
            cable(i)=0;
            delay(i)=276.4;
            delaycorr(i)=+5.1;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.4*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        %% Warning!!!! Inverted with 106
        case 105
            if nrun<=469100
                if gps
                    pos_ant(i,:)=[1742.1 -80.0 2656.0];
                else
                    pos_ant(i,:)=[1748.0  -81.0 2656.0];
                end
                podN(i) = 106;
                pos_pod(i,:)=[1800.0       0.0     	 2656.2];
                cable(i)=0;
                delay(i)=266.9;
                delaycorr(i)=0;
                if nruntrue>300000
                    delay(i) = delay(i)+ 0.4*PROPAGT;
                end
                sciID(i)=0;
                logID(i)=0; 
            else  % Fiber swapped with 106
                if gps
                    pos_ant(i,:)=[1742.4 69.8 2657.0];
                else
                    pos_ant(i,:)=[1748.0   68.0 2657.0];
                end
                podN(i) = 105;
                pos_pod(i,:)=[1980.0       0.0     	 2659.8];
                cable(i)=0;
                delay(i)=537.6;
                delaycorr(i)=-11;
                if nruntrue>300000
                    delay(i) = delay(i)+ 0.4*PROPAGT;
                end
                sciID(i)=0;
                logID(i)=0;                 
            end

        case 106
            if nrun<=469100
                if gps
                    pos_ant(i,:)=[1742.4 69.8 2657.0];
                else
                    pos_ant(i,:)=[1748.0   68.0 2657.0];
                end
                podN(i) = 105;
                pos_pod(i,:)=[1980.0       0.0     	 2659.8];
                cable(i)=0;
                delay(i)=537.6;
                delaycorr(i)=-11;
                if nruntrue>300000
                    delay(i) = delay(i)+ 0.4*PROPAGT;
                end
                sciID(i)=0;
                logID(i)=0; 
            else  % Fiber swapped with 105
                if gps
                    pos_ant(i,:)=[1742.1 -80.0 2656.0];
                else
                    pos_ant(i,:)=[1748.0  -81.0 2656.0];
                end
                podN(i) = 106;
                pos_pod(i,:)=[1800.0       0.0     	 2656.2];
                cable(i)=0;
                delay(i)=266.9;
                delaycorr(i)=0;
                if nruntrue>300000
                    delay(i) = delay(i)+ 0.4*PROPAGT;
                end
                sciID(i)=0;
                logID(i)=0; 
            end
            
        case 107
            if gps
                pos_ant(i,:)=[1892.7 -202.3 2659.0];
            else
                pos_ant(i,:)=[1899.0 -203.0 2659.0];
            end
            podN(i) = 107;
            pos_pod(i,:)=[2000.0       0.0     	 2660.4];
            cable(i)=0;
            delay(i)=547.0;
            delaycorr(i)=+6;
            sciID(i)=0;
            logID(i)=0; 
            
        case 108
            if gps
                pos_ant(i,:)=[1903.0 152.5 2659.0];
            else
                pos_ant(i,:)=[1908.5  151.0 2659.0];
            end
            podN(i) = 108;
            pos_pod(i,:)=[2100.0       0.0     	 2664.4];
            cable(i)=0;
            delay(i)=658.0;
            delaycorr(i)=+6;
            sciID(i)=0;
            logID(i)=0; 
            
        case 109
            if gps
                pos_ant(i,:)=[2226.7 -5.8 2668.0];
            else
                pos_ant(i,:)=[2233.0   -8.0 2668.0];
            end
            podN(i) = 109;
            pos_pod(i,:)=[2120.0       0.0     	 2665.2];
            cable(i)=0;
            delay(i)=584.2;
            delaycorr(i)=-27;
            sciID(i)=0;            
        
        case 110
            if gps
                pos_ant(i,:)=[2042.6 -81.1 2663.0];
            else
                pos_ant(i,:)=[2049.0  -82.8 2663.0];
            end
            podN(i) = 110;
            pos_pod(i,:)=[2140.0       0.0     	 2666.0];
            cable(i)=0;
            delay(i)=609.2;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 1.0*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 111
            if gps
                pos_ant(i,:)=[2042.6 68.3 2664.0];
            else
                pos_ant(i,:)=[2048.0   65.8 2664.0];
            end
            podN(i) = 111;
            pos_pod(i,:)=[2160.0       0.0     	 2666.8];
            cable(i)=0;
            delay(i)=647.0;
            delaycorr(i)=+22;
            sciID(i)=0;
            logID(i)=0; 
            
        case 112
            if gps
                pos_ant(i,:)=[2193.8 -222.3 2668.0];
            else
                pos_ant(i,:)=[2201.0 -224.4 2668.0];
            end
            podN(i) = 112;
            pos_pod(i,:)=[2180.0       0.0     	 2667.6];
            cable(i)=0;
            delay(i)=731.4;
            delaycorr(i)=+5;
            sciID(i)=0;
            logID(i)=0; 
            
        case 113
            if gps
                pos_ant(i,:)=[2197.6 188.7 2667.0];
            else
                pos_ant(i,:)=[2204.0  186.7 2667.0];
            end
            podN(i) = 113;
            pos_pod(i,:)=[2200.0       0.0     	 2668.1];
            cable(i)=0;
            delay(i)=728.3;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 114
            if gps
                pos_ant(i,:)=[2347.4 -81.1 2673.0];
            else
                pos_ant(i,:)=[2355.0  -82.0 2673.0];
            end
            podN(i) = 114;
            pos_pod(i,:)=[2300.0       0.0     	 2669.6];
            cable(i)=0;
            delay(i)=748.7;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.5*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 115
            if gps
                pos_ant(i,:)=[2348.9 68.5 2671.0];
            else
                pos_ant(i,:)=[2354.0   66.0 2671.0];
            end
            podN(i) = 115;
            pos_pod(i,:)=[2340.0       0.0     	 2670.2];
            cable(i)=0;
            delay(i)=785.7;
            delaycorr(i)=-6;
            sciID(i)=0;
            logID(i)=0; 
            
        case 116
            if gps
                pos_ant(i,:)=[2498.2 192.2 2674.0];
            else
                pos_ant(i,:)=[2505.0  191.0 2674.0];
            end
            podN(i) = 116;
            pos_pod(i,:)=[2380.0       0.0     	 2670.8];
            cable(i)=0;
            delay(i)=947.6;
            delaycorr(i)=-19;
            sciID(i)=0;
            logID(i)=0; 
            
        case 117
            if gps
                pos_ant(i,:)=[2495.7 -206.3 2684.0];
            else
                pos_ant(i,:)=[2504.0 -207.0 2684.0];
            end
            podN(i) = 117;
            pos_pod(i,:)=[2440.0       0.0     	 2672.7];
            cable(i)=0;
            delay(i)=977.1;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 118
            if gps
                pos_ant(i,:)=[2507.7 -6.4 2677.0];
            else
                pos_ant(i,:)=[2495.0   -8.0 2677.0];
            end
            podN(i) = 118;
            pos_pod(i,:)=[2700.0       0.0     	 2681.3];
            cable(i)=0;
            delay(i)=1206.5;
            delaycorr(i)=-11;
            sciID(i)=0;
            logID(i)=0; 
            
        case 119
            if gps
                pos_ant(i,:)=[2646.1 -82.0  2687.0];
            else
                pos_ant(i,:)=[2654.0 -82.3 2687.0];
            end
            podN(i) = 119;
            pos_pod(i,:)=[2820.0       0.0     	 2684.6];
            cable(i)=0;
            delay(i)=1320.0;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 120
            if gps
                pos_ant(i,:)=[2646.5 67.8 2680.0];
            else
                pos_ant(i,:)=[2655.0 67.0 2680.0];
            end
            podN(i) = 120;
            pos_pod(i,:)=[2840.0       0.0     	 2685.0];
            cable(i)=0;
            delay(i)=1358.2;
            delaycorr(i)=+16;
            sciID(i)=0;
            logID(i)=0; 
            
        case 121
            pos_ant(i,:)=[1319.0 -243.0 2646.8];
            podN(i) = 121;
            pos_pod(i,:)=[1380.0       0.0     	 2651.8];
            cable(i)=0;
            delay(i)=130.9;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 122
            if gps
                pos_ant(i,:)=[1299.7 187.1 2646.0];
            else
                pos_ant(i,:)=[1303.5  184.0 2646.0];
            end
            podN(i) = 122;
            pos_pod(i,:)=[1360.0       0.0     	 2651.4];
            cable(i)=0;
            delay(i)=182.2;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.5*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 123
            if gps
                pos_ant(i,:)=[1127.4 60.0 2645.0];
            else
                pos_ant(i,:)=[1128.0   57.0 2645.0];
            end
            podN(i) = 123;
            pos_pod(i,:)=[1240.0       0.0     	 2649.0];
            cable(i)=0;
            delay(i)=171.6;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 124
            if gps
                pos_ant(i,:)=[1133.2 -89.3 2645.0];
            else
                pos_ant(i,:)=[1135.0  -91.0 2645.0];
            end    
            podN(i) = 124;
            pos_pod(i,:)=[980.0       0.0     	 2643.8];
            cable(i)=0;
            delay(i)=470.7;
            delaycorr(i)=+3;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.5*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 125
            if gps
                pos_ant(i,:)=[983.6 170.4 2647.0];
            else
                pos_ant(i,:)=[983.0  166.0 2647.0];
            end
            podN(i) = 125;
            pos_pod(i,:)=[920.0       0.0     	 2642.6];
            cable(i)=0;
            delay(i)=535.5;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.6*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 126
            if gps
                pos_ant(i,:)=[989.3 -213.0 2647.0];
            else
                pos_ant(i,:)=[991.0 -216.0 2647.0];
            end
            podN(i) = 126;
            pos_pod(i,:)=[880.0       0.0     	 2641.8];
            cable(i)=0;
            delay(i)=612.1;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.55*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 127
            if gps
                pos_ant(i,:)=[848.5 -86.6 2644.0];
            else
                pos_ant(i,:)=[850.0  -90.0 2644.0];
            end
            podN(i) = 127;
            pos_pod(i,:)=[840.0       0.0     	 2641.0];
            cable(i)=0;
            delay(i)=539.0;
            delaycorr(i)=+5;
            sciID(i)=0;
            logID(i)=0; 
            
        case 128
            if gps
                pos_ant(i,:)=[848.4 62.0 2643.0];
            else
                pos_ant(i,:)=[849.0   58.0 2643.0];
            end
            podN(i) = 128;
            pos_pod(i,:)=[740.0       0.0	 2640.0];
            cable(i)=0;
            delay(i)=673.7;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 129
            if gps
                pos_ant(i,:)=[699.1 194.5 2643.0];
            else
                pos_ant(i,:)=[703.0  191.0 2643.0];
            end
            podN(i) = 129;
            pos_pod(i,:)=[720.0       0.0	 2640.0];
            cable(i)=0;
            delay(i)=740.0;
            delaycorr(i)=-14;
            sciID(i)=0;
            logID(i)=0; 
            
        case 130
            if gps
                pos_ant(i,:)=[699.4 -205.8 2643.0];
            else
                pos_ant(i,:)=[701.0 -210.0 2643.0];
            end
            podN(i) = 130;
            pos_pod(i,:)=[700.0       0.0	 2640.0 ];
            cable(i)=0;
            delay(i)=759.2;
            delaycorr(i)=+0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 131
            if gps
                pos_ant(i,:)=[549.3 66.3 2640.0];
            else
                pos_ant(i,:)=[549.0   64.0 2640.0];
            end
            podN(i) = 131;
            pos_pod(i,:)=[680.0       0.0	 2640.0];
            cable(i)=0;
            delay(i)=740.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 132
            if gps
                pos_ant(i,:)=[548.7 -79.9  2640.0];
            else
                pos_ant(i,:)=[549.0  -82.2 2640.0];
            end
            podN(i) = 132;
            pos_pod(i,:)=[660.0       0.0	 2640.0];
            cable(i)=0;
            delay(i)=746.4;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.5*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 133
            if gps
                pos_ant(i,:)=[479.9 165.6 2639.5];
            else
                pos_ant(i,:)=[480.0  163.0 2639.5];
            end
            podN(i) = 133;
            pos_pod(i,:)=[640.0       0.0	 2640.0];
            cable(i)=0;
            delay(i)=917.0;
            delaycorr(i)=0;
            if nruntrue>300000
                delay(i) = delay(i)+ 0.5*PROPAGT;
            end
            sciID(i)=0;
            logID(i)=0; 
            
        case 134
            if gps
                pos_ant(i,:)=[462.3 -130.6 2640.0];
            else
                pos_ant(i,:)=[462.0 -133.0 2640.0];
            end
            podN(i) = 134;
            pos_pod(i,:)=[540.0       0.0	 2639.2];
            cable(i)=0;
            delay(i)=874.3;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 135
            if gps
                pos_ant(i,:)=[369.6 80.6  2637.0];
            else
                pos_ant(i,:)=[369.0 78.5 2637.0];
            end
            podN(i) = 135;
            pos_pod(i,:)=[520.0       0.0	 2638.9];
            cable(i)=0;
            delay(i)=909.1;
            delaycorr(i)=0;
            sciID(i)=0;            
            logID(i)=0; 
            
        case 136
            pos_ant(i,:)=[304.0 -112.0 2636.0];
            podN(i) = 136;
            pos_pod(i,:)=[340.0       0.0	 2636.2];
            cable(i)=0;
            delay(i)=1101.0;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 137
            if gps
                pos_ant(i,:)=[232.0 26.6  2633.0];
            else
                pos_ant(i,:)=[231.0   23.0 2633.0];
            end
            podN(i) = 137;
            pos_pod(i,:)=[260.0       0.0	 2635.0];
            cable(i)=0;
            delay(i)=1064.9;
            delaycorr(i)= -7;
            sciID(i)=0;
            logID(i)=0; 
            
        case 138
            if gps
                pos_ant(i,:)=[224.6 -45.3 2634.0];
            else
                pos_ant(i,:)=[224.0  -49.0 2634.0];
            end
            podN(i) = 138;
            pos_pod(i,:)=[240.0       0.0	 2634.7];
            cable(i)=0;
            delay(i)=1098.6;
            delaycorr(i)=-8;
            sciID(i)=0;
            logID(i)=0; 
            
        case 140
            if gps
                pos_ant(i,:)=[100.6 -104.6 2635.0];
            else
                pos_ant(i,:)=[100.0 -108.0 2635.0];
            end
            podN(i) = 140;
            pos_pod(i,:)=[100.0       0.0     	 2633.1];
            cable(i)=0;
            delay(i)=1257.9;
            delaycorr(i)=+4;
            sciID(i)=0;
            logID(i)=0; 
            
        case 148
            %pos_ant(i,:)=[61.4 566.7 2645.1];
            pos_ant(i,:)=[49.9, 559.5, 2640.5];
            podN(i) = 148;
            pos_pod(i,:)=[0.0	   527.0     	 2638.7];
            cable(i)=0;
            delay(i)=1893.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 149
            %pos_ant(i,:)=[-18.1 532.1 2644.6];
            pos_ant(i,:)=[-19.9, 530.3, 2640.2];
            podN(i) = 149;
            pos_pod(i,:)=[0.0	   497.7     	 2639.0];
            cable(i)=0;
            delay(i)=1841.6;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 150
            %pos_ant(i,:)=[100.3 467.7 2641.9];
            pos_ant(i,:)=[90.8, 460.2, 2638.3];
            podN(i) = 150;
            pos_pod(i,:)=[0.0	   468.3     	 2639.3];
            cable(i)=0;
            delay(i)=1828.8;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 151
            %pos_ant(i,:)=[35.0 457.0 2642.0];
            pos_ant(i,:)=[23.7, 449.4, 2639.4]; % Jan 2012
            podN(i) = 151;
            pos_pod(i,:)=[0.0	   438.9     	 2639.6 ];
            cable(i)=0;
            delay(i)=1756.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 152
            %pos_ant(i,:)=[-22.7 394.5 2641.9];
            pos_ant(i,:)=[-33.4, 385.6, 2637.4]; % Jan 2012
            podN(i) = 152;
            pos_pod(i,:)=[0.0	   409.6     	 2639.9];
            cable(i)=0;
            delay(i)=1728.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 153
            %pos_ant(i,:)=[92.0 349.5 2638.9];
            pos_ant(i,:)=[82.8, 342.0, 2635.7]; % Jan 2012
            podN(i) = 153;
            pos_pod(i,:)=[0.0	   380.2     	 2639.3];
            cable(i)=0;
            delay(i)=1740.2;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 154
            %pos_ant(i,:)=[-24.8 263.9 2634.9];
            pos_ant(i,:)=[-33.6, 256.6, 2633.1]; % Jan 2012
            podN(i) = 154;
            pos_pod(i,:)=[0.0	   233.4     	 2632.7];
            cable(i)=0;
            delay(i)=1553.3;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 155
            %pos_ant(i,:)=[66.7 196.6 2635.3]; %Meng
            pos_ant(i,:)=[58.1, 190.1, 2632.0]; % Jan 2012
            podN(i) = 155;
            pos_pod(i,:)=[0.0	   204.0     	 2631.4];
            cable(i)=0;
            delay(i)=1535.0;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 156
            %pos_ant(i,:)=[-11.3 11.3 2631.5];  %21CMA antenna (removed)
            pos_ant(i,:)=[-20.6 -1 2632.7];
            podN(i) = 156;
            pos_pod(i,:)=[0.0	  -60.3     	 2633.4];
            cable(i)=0;
            delay(i)=1400;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=0; 
            
        case 157
            %pos_ant(i,:)=[21.2 -95.6 2635.8];
            pos_ant(i,:)=[22.5, -100.7, 2634.8];
            podN(i) = 157;
            pos_pod(i,:)=[0.0	  -177.8     	 2636.3];
            cable(i)=0;
            delay(i)=1515.5;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1; 
            
        case 158
            %pos_ant(i,:)=[81.1 -191.7 2637.7];
            %pos_ant(i,:)=[72 -197 2637.7];
            pos_ant(i,:)=[72 -199 2637.7];
            podN(i) = 158;
            pos_pod(i,:)=[0.0	  -207.1     	 2637.1];
            cable(i)=0;
            delay(i)=1563.4;
            delaycorr(i)=0;
            sciID(i)=0;
            logID(i)=1;     
        
        %% Scintillators
        case 139  
            pos_ant(i,:)=[129.0  -34.4 2634.1];
            podN(i) = 139;
            pos_pod(i,:)=[140.0       0.0	 2633.5];
            cable(i)=0;
            delay(i)=1194; %1233.2;
            delaycorr(i)=0;
            sciID(i)=1;
            logID(i)=0;            
            
        case 175
            pos_ant(i,:)=[4.5     191.4  2631.4];
            podN(i) = 155;
            pos_pod(i,:)=[0.0	   204.0     	 2631.4];
            cable(i)=0;
            delay(i)=1505.8; %1535;
            delaycorr(i)=0;
            sciID(i)=1;            
            logID(i)=0; 
            
        case 176
            pos_ant(i,:)=[0.0     -37.5  2633.4];
            podN(i) = 156;
            pos_pod(i,:)=[0.0	  -60.3     	 2633.4];
            cable(i)=0;
            delay(i)=1346.4; %1400;
            delaycorr(i)=0;
            sciID(i)=1;
            logID(i)=0; 
            
        otherwise
            display(sprintf('SetupCharacteristics error : Unknow antenna %3d',ant(j)))
            break;
    end;
    
end;
