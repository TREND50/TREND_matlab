function [Struct]=ComputeBarycenter(Struct,slice)

SharedGlobals;

nrun = Struct.Setup.Run;
ncoinc=Struct.Setup.TotalCoinc;
isScint = [Struct.Setup.Det.isScint];
Lant = Struct.Coinc.MultAnt;
X=[Struct.Setup.Det.X];
Y=[Struct.Setup.Det.Y];
Z=[Struct.Setup.Det.Z];
status = Struct.Coinc.Det.Status;
trig=[status>=0];
if RawFilterFlag == 1
    in = [status>=0 & status~=2];
else
    in = trig;
end
in = [Struct.Coinc.Det.Evt>0];
dm = zeros(ncoinc,1);

for i=1:ncoinc
    ini=in(i,isScint==0);
    if sum(ini)<4
        continue
    end
    % Barycenter position
    xb = mean(X(ini));
    yb = mean(Y(ini));
    zb = mean(Z(ini));
    %
    % Mean Distance to barycenter
    if Lant(i) == sum(ini)
        dis = sqrt( (X(ini)-xb*ones(1,Lant(i))).^2 + (Y(ini)-yb*ones(1,Lant(i))).^2 + (Z(ini)-zb*ones(1,Lant(i))).^2 );
        dm(i) = mean(dis);
    else
        disp(sprintf('Warning! Coinc %d: Mult = %d while nb events = %d',i,Lant(i),sum(ini)));
        dm(i) = 1e8;
    end
    %pause
end;	

Struct.Coinc.DistBary = dm;        
        
filename =['D:/dst/cleanfilt/' sprintf(dst_filename,nrun,slice)];
save(filename,'Struct');
display(sprintf('DST %s now saved to file.',filename));
