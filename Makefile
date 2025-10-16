# Executable.
NAME = matins

# Paths
BASE = $(CURDIR)
EXEC = $(BASE)/$(NAME)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
SRCDi = $(BASE)/src/initial
SRCDm = $(BASE)/src/magnetic
SRCDp = $(BASE)/src/microphysics
SRCDo = $(BASE)/src/output
SRCDt = $(BASE)/src/thermal

# Compiler, choose (INTEL_MKL, INTEL_openblas, GNU).
COMPILATION=GNU
ifeq ($(COMPILATION),INTEL_openblas)
  FCOMP = ifort
  FFLAGS = -O3 -xHost -parallel -ipo -fno-math-errno -fp-model precise -fp-model source -ftz -heap-arrays 1024 -qopenmp -check bounds -g -traceback
  FLINK = -llapack -liomp5 -lpthread -lm -ldl
endif
ifeq ($(COMPILATION),INTEL_MKL)
  FCOMP = ifort
  FFLAGS = -O3 -xHost -parallel -ipo -fno-math-errno -fp-model precise -fp-model source -ftz -heap-arrays 2048 -qopenmp -qmkl=parallel -check bounds -g -traceback
  FLINK = -liomp5 -lpthread -lm -ldl
endif
ifeq ($(COMPILATION),GNU)
  FCOMP = gfortran
  FFLAGS = -ffpe-trap=invalid,zero -fopenmp -O3 -fbounds-check # -pg
  FLINK = -llapack
endif

# Search paths.
VPATH = $(OBJD):$(SRCD):$(SRCDi):$(SRCDm):$(SRCDp):$(SRCDo):$(SRCDt)

# Object list.
OBJS = $(OBJD)/main.o $(OBJD)/constants.o $(OBJD)/grid.o $(OBJD)/initial.o $(OBJD)/initial_magnetic.o $(OBJD)/input.o $(OBJD)/ns_structure.o $(OBJD)/odeint.o $(OBJD)/bevol.o $(OBJD)/magnetic_analysis.o $(OBJD)/cond_crust_potekhin19.o $(OBJD)/cond_inner_crust_ei.o $(OBJD)/cond_interface.o $(OBJD)/cond_outer_crust_ei.o $(OBJD)/cond_phonons.o $(OBJD)/effective_mass.o $(OBJD)/emissivity.o $(OBJD)/eos14.o $(OBJD)/eosmag14.o $(OBJD)/get_Zimp.o $(OBJD)/heat_capacity.o $(OBJD)/simple_microphysics.o $(OBJD)/superfluid.o $(OBJD)/output.o $(OBJD)/utils.o $(OBJD)/envelope.o $(OBJD)/tevol.o $(OBJD)/timestep_cooling.o

# Linking the object files to build an executable.
$(NAME): $(OBJS)
	@echo Building executable...
	$(FCOMP) -o $@ $(OBJS) $(FFLAGS) $(PRECISION) $(FLINK)
	@echo All done...
	@echo Working directory: $(CURDIR)
# If the working directory is not correct, the paths must be edited.

# Compiling objects.

# main
$(OBJD)/main.o: $(SRCD)/main.f90 $(OBJD)/output.o $(OBJD)/tevol.o $(OBJD)/bevol.o $(OBJD)/cond_interface.o $(OBJD)/simple_microphysics.o $(OBJD)/grid.o $(OBJD)/input.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $<

# initial
$(OBJD)/constants.o: $(SRCDi)/constants.f90
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/input.o: $(SRCDi)/input.f90
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/grid.o: $(SRCDi)/grid.f90 $(OBJD)/constants.o $(OBJD)/input.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/ns_structure.o: $(SRCDi)/ns_structure.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/initial_magnetic.o: $(SRCDi)/initial_magnetic.f90 $(OBJD)/bevol.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/initial.o: $(SRCDi)/initial.f90 $(OBJD)/grid.o $(OBJD)/constants.o $(OBJD)/input.o $(OBJD)/initial_magnetic.o $(OBJD)/bevol.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

# output
$(OBJD)/utils.o: $(SRCDo)/utils.f90
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/output.o: $(SRCDo)/output.f90 $(OBJD)/constants.o $(OBJD)/grid.o $(OBJD)/utils.o $(OBJD)/magnetic_analysis.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

# thermal
$(OBJD)/envelope.o: $(SRCDt)/envelope.f90 $(OBJD)/grid.o $(OBJD)/constants.o $(OBJD)/input.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/tevol.o: $(SRCDt)/tevol.f90 $(OBJD)/grid.o $(OBJD)/constants.o $(OBJD)/simple_microphysics.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/timestep_cooling.o: $(SRCDt)/timestep_cooling.f90 $(OBJD)/input.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

# magnetic
$(OBJD)/bevol.o: $(SRCDm)/bevol.f90 $(OBJD)/grid.o $(OBJD)/input.o $(OBJD)/magnetic_analysis.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/magnetic_analysis.o: $(SRCDm)/magnetic_analysis.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

# mycrophysics
$(OBJD)/cond_interface.o: $(SRCDp)/cond_interface.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/effective_mass.o: $(SRCDp)/effective_mass.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/emissivity.o: $(SRCDp)/emissivity.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/get_Zimp.o: $(SRCDp)/get_Zimp.f90 $(OBJD)/input.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/heat_capacity.o: $(SRCDp)/heat_capacity.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/simple_microphysics.o: $(SRCDp)/simple_microphysics.f90 $(OBJD)/constants.o $(OBJD)/grid.o
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/superfluid.o: $(SRCDp)/superfluid.f90
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $<

$(OBJD)/%.o: $(SRCDi)/%.f
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 

$(OBJD)/%.o: $(SRCDp)/%.f
	@echo Building $*.o...
	cd $(OBJD); $(FCOMP) -c $(FFLAGS) $(PRECISION) $< 


# Other rules.
clean:
	rm -f -r $(EXEC) obj/*.o obj/*.mod ./out/1D/*.yg ./out/1D/*.d ./out/2D/*.dat ./out/3D/*.vtu ./out/3D/*.pvd ./out/energy/*.dat ./out/energy/*.yg ./out/energy/*.d

outremove:
	rm -f -r ./out/1D/*.yg ./out/1D/*.d ./out/2D/*.dat ./out/3D/*.vtu ./out/3D/*.pvd ./out/energy/*.dat ./out/energy/*.yg ./out/energy/*.d

out_temp:
	mkdir OUT_TEMP
	cp in/input.dat OUT_TEMP/
	mkdir OUT_TEMP/1D
	mv out/1D/*.yg OUT_TEMP/1D/
	mv out/1D/*.d OUT_TEMP/1D/
	mkdir OUT_TEMP/2D
	mv out/2D/*.dat OUT_TEMP/2D/
	mkdir OUT_TEMP/3D
	mv out/3D/*.vtu OUT_TEMP/3D/
	mv out/3D/*.pvd OUT_TEMP/3D/
	mkdir OUT_TEMP/energy
	mv out/energy/*.dat OUT_TEMP/energy/
	mv out/energy/*.yg OUT_TEMP/energy/
	cp -r out/python OUT_TEMP/

copy:
	rm -rf ../code_bckp
	mkdir ../code_bckp
	cp makefile ../code_bckp/.
	cp -r out ../code_bckp/.
	cp -r obj ../code_bckp/.
	cp -r in ../code_bckp/.
	cp -r src ../code_bckp/.

