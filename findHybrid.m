function [] = findHybrid()

SharedGlobals

nrun=2550:2710;

for i=1:length(nrun)
    irun = nrun(i);
    for ifile = 1 :3
        dstname = [DST_PATH sprintf(dst_filename,irun,ifile)];
        if ~exist(dstname)
            %disp(sprintf('dst %s not found',dstname))
            continue
        end
        dst = load(dstname);
        format long
        (max(dst.Struct.Setup.InfosRun.TimeStop)- min(dst.Struct.Setup.InfosRun.TimeStart))/60;
        dur = (max(dst.Struct.Setup.RunTimeStop)-min(dst.Struct.Setup.RunTimeStart))/60;
        
        
        ls = dst.Struct.Coinc.MultSci;
        la = dst.Struct.Coinc.MultAnt;
        idcoinc = dst.Struct.Coinc.IdCoinc;
        det = dst.Struct.Coinc.Det.Id;
        evt = dst.Struct.Coinc.Det.Evt;
        tag = dst.Struct.Coinc.Det.Tag;
        hyb = find ( ls>1 & ls+la>3);
        nhyb = size(hyb,1);
        type = [dst.Struct.Setup.Det.isScint];

        if nhyb==0
            disp(sprintf('No hybrid in run %d',irun))
        else
            for j = 1 :nhyb
                disp(sprintf('R%d coinc %d: %d scints + %d antennas',irun,idcoinc(hyb(j)),ls(hyb(j)),la(hyb(j))))
                disp 'Scints:'
                Id = det(hyb(j),tag(hyb(j),:)==1 & type==1)
                Evt = evt(hyb(j),tag(hyb(j),:)==1 & type==1)
                disp 'Antennas:'
                Id = det(hyb(j),tag(hyb(j),:)==1 & type==0)
                Evt = evt(hyb(j),tag(hyb(j),:)==1 & type==0)
                
            end
        end
        clear dst
    end
end
