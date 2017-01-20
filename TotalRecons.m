function [baddst] = TotalRecons(runstart,runstop)
% Compute & stores zenith & azim angle distributions
% for ALL reconstructed events.
% OMH 29/11/2012


SharedGlobals;
ThetaSTot = zeros(1,100);
PhiSTot = zeros(1,360);
ThetaPTot = zeros(1,100);
PhiPTot = zeros(1,360);
Chi2PTot = zeros(1,1000);
Chi2STot = zeros(1,1000);
RadTot = zeros(1,100);
RampTot = zeros(1,200);
MultTot = zeros(1,50);

nbins = 250;
XTot = zeros(1,nbins);
YTot = zeros(1,nbins);
ZTot = zeros(1,nbins);
XYTot = zeros(nbins,nbins);
Runs = [];
Trigs = [];
CoincRaw = [];
CoincConsCoinc = [];
CoincBadPulses = [];
Mult5 = [];
Valid = [];
baddst = [];

for i=runstart:runstop
    %% Get number of ub dsts
    stopflag=0;
    nbiter=1;
    while stopflag==0
        if LIGHT 
          filename = [DST_PATH sprintf('dst%d_%d_light.mat',i,nbiter)];
        else
          filename = [DST_PATH sprintf('dst%d_%d.mat',i,nbiter)];
        end
        fd=fopen(filename);
        if fd~=-1
            nbiter=nbiter+1;
            fclose(fd);
        else
            stopflag=1;
        end;
    end;
    nbiter=nbiter-1;
    
    if nbiter==0
        display(sprintf('No dst found for run %d',i))
        continue
    end;
    
    %% Loop on meta dsts
    for j=1:nbiter       
        if LIGHT 
            filename = [DST_PATH sprintf('dst%d_%d_light.mat',i,j)];
        else
             filename = [DST_PATH sprintf('dst%d_%d.mat',i,j)];
        end
        display(sprintf('Loading %s...',filename))
        d = load(filename);
        display 'Done.'
        
       if ~isfield(d.Struct.Setup.InfosRun,'MetaTotalCoinc')
            disp(sprintf('No field Candidates in %s!',filename));
            baddst(end+1) = i;
            continue
        end

        %% Load dst
        NEvt = d.Struct.Setup.TotalEvt;
        NCoincRaw = sum(d.Struct.Setup.InfosRun.MetaTotalCoinc);
        NCoinc_ConsCoinc = sum(d.Struct.Setup.InfosRun.TotalCoincFiltered1);
        NCoinc_BadPulses = d.Struct.Setup.TotalCoinc;
        Mult = d.Struct.Coinc.MultAnt';
        if CORDELAYS
            PlanStruct = d.Struct.Coinc.DelayCorrRecons.PlanRecons.Radio;
            SphStruct = d.Struct.Coinc.DelayCorrRecons.SphRecons;
        else
            PlanStruct = d.Struct.Coinc.PlanRecons.Radio;
            SphStruct = d.Struct.Coinc.SphRecons;
        end
        Amp = d.Struct.Coinc.Det.AmpMax;
        maxAmp = max(Amp,[],2);
        Amp(Amp==0) = 1e8;
        minAmp = min(Amp,[],2);
        ramp = maxAmp./minAmp;
        ThetaP = PlanStruct.Theta;
        PhiP = PlanStruct.Phi;
        Slopep = PlanStruct.SlopeDelay;
        Chi2p = PlanStruct.Chi2Delay;
        x0 = SphStruct.X0;
        y0 = SphStruct.Y0;
        z0 = SphStruct.Z0;
        ThetaS = SphStruct.Theta;
        PhiS = SphStruct.Phi;
        Radius = SphStruct.minDistSource;
        Slopes = SphStruct.SlopeDelay;
        Chi2s = SphStruct.Chi2Delay;
        
        %% Cuts
        all = find(Mult>=5);
        okfar = find( abs(Slopes-1)<0.1 & Chi2s<50 & abs(Slopep-1)<0.1 & Chi2p<50 & Mult>=5 & Radius>=500);
        okvclose = find( abs(Slopes-1)<0.1 & Chi2s<50 & Mult>=5 & Radius<=500);
        okclose = find( abs(Slopes-1)<0.1 & Chi2s<50 & Mult>=5 & Radius<=500);
        %
        valid = [okvclose okfar];
        %
        disp(sprintf('Fraction of coincs passing cut 1 (ConsCoinc) = %d/%d (%3.1f pc).',NCoinc_ConsCoinc,NCoincRaw,NCoinc_ConsCoinc/NCoincRaw*100))
        disp(sprintf('Fraction of coincs passing cut 2 (BadPulses) = %d/%d (%3.1f pc).',NCoinc_BadPulses,NCoinc_ConsCoinc,NCoinc_BadPulses/NCoinc_ConsCoinc*100))
        disp(sprintf('Fraction of coincs reconstructed succesfully with L>4 = %d/%d (%3.1f pc).',length(valid),length(all),length(valid)/length(all)*100))
        disp(sprintf('Fraction of distant events = %d/%d (%3.1f pc)',length(okfar),length(valid),length(okfar)/length(valid)*100))
        
        %% Fill in histos
        %
        %Chi2
        for k = 1:1000
          %lChi2s = log10(Chi2s);  
          sel = find(Chi2s(all)>=k-1 & Chi2s(all)<k);
          Chi2STot(k) = Chi2STot(k) + length(sel);
        end
        for k = 1:1000
          %lChi2p = log10(Chi2p);  
          sel = find(Chi2p(all)>=k-1 & Chi2p(all)<k);
          Chi2PTot(k) = Chi2PTot(k) + length(sel);
        end
        for k = 1:100
          lradius = log10(Radius);
          sel = find(lradius(valid)>=(k-1)/10 & lradius(valid)<k/10);
          RadTot(k) = RadTot(k) + length(sel);
        end
        for k = 1:50
          sel = find(Mult(valid)>=k-1 & Mult(valid)<k);
          MultTot(k) = MultTot(k) + length(sel);
        end
        for k = 1:200
          sel = find(ramp(valid)>=(k-1)/10 & ramp(valid)<k/10);
          RampTot(k) = RampTot(k) + length(sel);
        end
        
        % 
        % Distant sources
        for k = 1:100
          sel = find(ThetaS(okfar)>=k-1 & ThetaS(okfar)<k);
          ThetaSTot(k) = ThetaSTot(k) + length(sel);
          sel = find(ThetaP(okfar)>=k-1 & ThetaP(okfar)<k);
          ThetaPTot(k) = ThetaPTot(k) + length(sel);
        end
        for k = 1:360
          sel = find(PhiS(okfar)>=k-1 & PhiS(okfar)<k);
          PhiSTot(k) = PhiSTot(k) + length(sel);
          sel = find(PhiP(okfar)>=k-1 & PhiP(okfar)<k);
          PhiPTot(k) = PhiPTot(k) + length(sel);
        end   
        
        % Close sources
        xmin = -1000;
        ymin = -1500;
        zmin = 2000;
        step = 20;
        for k = 1:nbins
          xa = xmin+(k-1)*step;
          xb = xmin+k*step;
          sel = find(x0(okclose)>=xa & x0(okclose)<=xb);
          selxy = find(x0(valid)>=xa & x0(valid)<=xb);
          XTot(k) = XTot(k) + length(sel);
          for kk = 1:nbins
            ya = ymin+(kk-1)*step;
            yb = ymin+kk*step;
            sel2 = find(y0(valid)>=ya & y0(valid)<=yb);
            XYTot(k,kk) = XYTot(k,kk)+length(intersect(selxy,sel2));
          end
          ya = ymin+(k-1)*step;
          yb = ymin+k*step;
          sel = find(y0(okclose)>=ya & y0(okclose)<=yb);
          YTot(k) = YTot(k) + length(sel);
          za = zmin+(k-1)*step;
          zb = zmin+k*step;
          sel = find(z0(okclose)>=za & z0(okclose)<=zb);
          ZTot(k) = ZTot(k) + length(sel);
        end
        
        Runs(end+1) = i;
        Trigs(end+1) = NEvt;
        CoincRaw(end+1) = NCoincRaw;
        CoincConsCoinc(end+1) = NCoinc_ConsCoinc;
        CoincBadPulses(end+1) = NCoinc_BadPulses;
        Mult5(end+1) = length(all);
        Valid(end+1) = length(valid);

     end
end

%% Save to file
filename = sprintf('TotalRecons_R%dR%d_102014.mat',runstart,runstop)
save(filename,'Runs','Trigs','CoincRaw','CoincConsCoinc','CoincBadPulses','Mult5','Valid','Chi2STot','Chi2PTot','ThetaPTot','PhiPTot','ThetaSTot','PhiSTot','XTot','YTot','ZTot','XYTot','RadTot','RampTot','MultTot');

%% Plots 
% figure(1)
% plot(ThetaSTot,'sk')
% grid on
% figure(2)
% plot(PhiSTot,'sk')
% grid on
