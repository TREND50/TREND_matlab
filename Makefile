target := DstBuilder
binary := $(target)_centos7
#mcc    := /usr/local/matlab2012b/bin/mcc
#mcc    := /usr/local/$(MATLAB_VER)/bin/mcc
mcc := /usr/local/matlab/R2016b/bin/mcc
$(THRONG_DIR)/soft/bin/$(binary): matlab_compil/$(target).m
	rm -rf matlab_compil
	mkdir matlab_compil
	cp *.m matlab_compil/.
	cp ../matlab_tools/*.m matlab_compil/.
	$(mcc) -m -v -o $(binary) $<
	@echo "Copying binary"
	ln -s $(THRONG_DIR)/soft/ana/TREND_matlab/$(binary) ../../bin/.
	@echo "Cleaning C files"
	rm -rf $(target).prj $(target)_main.c $(target)_mcc_component_data.c mccExcludedFiles.log readme.txt run_$(target).sh
	@echo "All done"

clean:
	@echo "Cleaning C files"
	rm -rf $(target).prj $(target)_main.c $(target)_mcc_component_data.c mccExcludedFiles.log readme.txt run_$(target).sh 
	@echo "Cleaning binary"
	rm -rf ../../bin/$(binary) ./$(binary)
	@echo "All clean"
