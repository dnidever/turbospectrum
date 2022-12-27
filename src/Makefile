ifneq (,$(shell which ifort))
type = ifort
else ifneq (,$(shell which gfortran))
type = gfortran
else
type = error
endif

all: $(value type)
ifort:
	#-rm DATA
	#ln -fs DATA DATA
	(cd exec; make)
	cp -f exec/*_lu ../bin/
gfortran:
	#-rm DATA
	#ln -fs DATA DATA
	(cd exec-gf; make)
	cp -f exec-gf/*_lu ../bin/
clean:
	rm -f *.o
	rm -f exec/*.o exec/*_lu
	rm -f exec-gf/*.o exec-gf/*_lu
error:
	@echo "ERROR: No ifort or gfortran to compile the code"
