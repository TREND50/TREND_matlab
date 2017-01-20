function [Struct]=Dist2Source(Struct)

ncoinc=Struct.Setup.TotalCoinc;
Tag=Struct.Coinc.Det.Tag;
Detectors=[Struct.Setup.Det.Name];
indant=find([Struct.Setup.Det.isScint]==0);
X=[Struct.Setup.Det.X];
Y=[Struct.Setup.Det.Y];
Z=[Struct.Setup.Det.Z];
DistSource=-1*ones(size(Tag,1),size(Tag,2));
minDistSource = zeros(ncoinc,1);

for i=1:ncoinc
    indanttag=indant(Tag(i,indant)==1);

	if length(indanttag)>3
        if ~isfield(Struct.Coinc,'DelayCorrRecons')
            X0=Struct.Coinc.SphRecons.X0;
            Y0=Struct.Coinc.SphRecons.Y0;
            Z0=Struct.Coinc.SphRecons.Z0;
         else
             X0=Struct.Coinc.DelayCorrRecons.SphRecons.X0;
             Y0=Struct.Coinc.DelayCorrRecons.SphRecons.Y0;
             Z0=Struct.Coinc.DelayCorrRecons.SphRecons.Z0;
        end;
    	for j=1:length(indanttag)
        	% distance antenna source
        	DistSource(i,indanttag(j))=sqrt((X(indanttag(j))-X0(i))^2+(Y(indanttag(j))-Y0(i))^2+(Z(indanttag(j))-Z0(i))^2);
    	end;
    	minDistSource(i) = min(DistSource(i,indanttag));
	end;
end;	

if ~isfield(Struct.Coinc,'DelayCorrRecons')
    Struct.Coinc.SphRecons.DistSource=DistSource;
    Struct.Coinc.SphRecons.minDistSource=minDistSource;
else
    Struct.Coinc.DelayCorrRecons.SphRecons.DistSource=DistSource;
    Struct.Coinc.DelayCorrRecons.SphRecons.minDistSource=minDistSource;
end;

        
        
