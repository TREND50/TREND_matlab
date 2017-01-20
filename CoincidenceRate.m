function [GlobalCoincRate,DetCoincRate]=CoincidenceRate(Struct)

if Struct.Setup.TotalCoinc==0
   GlobalCoincRate=0;
   DetCoincRate=0;
   return;
end;

timescale=[Struct.Setup.InfosRun.TrigTime];
timeevt=max(Struct.Coinc.Det.UnixTime,[],2);
timescale = timescale + timeevt(1)-timescale(1);
det=[Struct.Setup.Det.Name];
tag=Struct.Coinc.Det.Tag;

GlobalCoincRate=zeros(size([Struct.Setup.InfosRun.TrigTime],1),1);
DetCoincRate=zeros(size([Struct.Setup.InfosRun.TrigRate],1),size([Struct.Setup.InfosRun.TrigRate],2));

% Global coincidence rate
for i=1:length(timescale)-1  
    ind=find(timeevt>=timescale(i) & timeevt<timescale(i+1));
    GlobalCoincRate(i)=length(ind)/(timescale(i+1)-timescale(i));
    clear ind
end;
GlobalCoincRate(end)=length(find(timeevt>=timescale(end)))/(Struct.Setup.RunTimeStop-timescale(end));

% Antenna coincidence rate
for i=1:length(det)
    indtag=find(tag(:,i)==1);
    for j=1:length(timescale)-1  
        ind=find(timeevt(indtag)>=timescale(j) & timeevt(indtag)<timescale(j+1));
        DetCoincRate(j,i)=length(ind)/(timescale(j+1)-timescale(j));
        clear ind
    end;
    DetCoincRate(end,i)=length(find(timeevt(indtag)>=timescale(end)))/(Struct.Setup.RunTimeStop-timescale(end));
    clear indtag
end;
