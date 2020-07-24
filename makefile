# makefile for themis program

# program name
PROGRAM = themis


# version control 
# Adding revision informatio to themis

#TODO:
#REV = $(shell git rev-parse --short HEAD)
# an error occurs

FC = gfortran

#FFLAGS = -g -Wall -O3 -ffpe-trap=invalid,zero,overflow -fcheck=bounds -Wno-compare-real -Wno-conversion -fbacktrace -fcheck=all -Wextra 
FFLAGS = -g -Wall -O3 -fcheck=bounds -Wno-compare-reals -Wno-conversion -fbacktrace -fcheck=all -Wextra -no-pie -ffpe-summary=none

# source directory
SRCDIR = src

# source files and objects
SRC = $(wildcard $(SRCDIR)/themis.f90)

# objects created
OBJS = mod_error_handling.o \
	mod_constants.o \
	mod_info.o \
	mod_cmd_line.o \
	mod_inquire.o \
	mod_input_read.o \
	mod_read_molecules.o \
	mod_spherical_grids.o \
	mod_grids.o \
	mod_pot_ljc.o \
	mod_pot_bhc.o \
	mod_pot_ljc_pair.o \
	xdr.o \
	libxdrfile.a \
	mod_loops.o \
	mod_write_vmd.o \
	mod_search_structures.o \
	mod_resume.o \
	mod_deallocate_all.o

# modules created
MODS =	mod_error_handling.mod \
	mod_constants.mod \
	mod_info.mod \
	mod_constants.mod \
	mod_cmd_line.mod \
	mod_inquire.mod \
	mod_input_read.mod \
	mod_read_molecules.mod \
	mod_spherical_grids.mod \
	mod_grids.mod \
	mod_pot_ljc.mod \
	mod_pot_bhc.mod \
	mod_pot_ljc_pair.mod \
	xdr.mod \
	mod_loops.mod \
	mod_write_vmd.mod \
	mod_search_structures.mod \
	mod_resume.mod \
	mod_deallocate_all.mod

all: $(OBJS)
	@printf "\n"
	@printf "\n"
	@printf "  GENERATING EXECUTABLE...\n"
	@printf "\n"
	@$(FC) $(SRC) $(FFLAGS) $(OBJS) -o $(PROGRAM) 
	@printf "   ** Compiling themis\n"
	@printf "\n"
	@printf "  CLEANING MODULES...\n"
	@printf "\n"
	@printf "   ** Deleting *.mod files\n"
	@rm -rf $(MODS) 
	@printf "\n"
	@printf "  CLEANING OBJECTS...\n"
	@printf "\n"
	@printf "   ** Deleting *.o files\n"
	@rm -rf $(OBJS)
	@printf "\n"

mod_error_handling.o: $(SRCDIR)/mod_error_handling.f90
	@printf "\n"
	@printf "  CREATING OBJECTS...\n"
	@printf "\n"
	@$(FC) -c $(SRCDIR)/mod_error_handling.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"       	

mod_constants.o: $(SRCDIR)/mod_constants.f90
	@$(FC) -c $(SRCDIR)/mod_constants.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_info.o: $(SRCDIR)/mod_info.f90
	@$(FC) -c $(SRCDIR)/mod_info.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"       	

mod_cmd_line.o: $(SRCDIR)/mod_cmd_line.f90
	@$(FC) -c $(SRCDIR)/mod_cmd_line.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_inquire.o: $(SRCDIR)/mod_inquire.f90
	@$(FC) -c $(SRCDIR)/mod_inquire.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_input_read.o: $(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_constants.f90
	@$(FC) -c $(SRCDIR)/mod_input_read.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_read_molecules.o: $(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_cmd_line.f90
	@$(FC) -c $(SRCDIR)/mod_read_molecules.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_spherical_grids.o: $(SRCDIR)/mod_spherical_grids.f90 \
	$(SRCDIR)/mod_constants.f90
	@$(FC) -c $(SRCDIR)/mod_spherical_grids.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_grids.o: $(SRCDIR)/mod_grids.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_cmd_line.f90 \
	$(SRCDIR)/mod_read_molecules.f90
	@$(FC) -c $(SRCDIR)/mod_grids.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_pot_ljc.o: $(SRCDIR)/mod_pot_ljc.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90
	@$(FC) -c $(SRCDIR)/mod_pot_ljc.f90 $(FFLAGS) 
	@printf "   ** Compiling $@\n"

mod_pot_bhc.o: $(SRCDIR)/mod_pot_bhc.f90 \
	$(SRCDIR)/mod_input_read.f90 \
  $(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90
	@$(FC) -c $(SRCDIR)/mod_pot_bhc.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_pot_ljc_pair.o: $(SRCDIR)/mod_pot_ljc_pair.f90 \
	$(SRCDIR)/mod_input_read.f90 \
  $(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90
	@$(FC) -c $(SRCDIR)/mod_pot_ljc_pair.f90 $(FFLAGS) 
	@printf "   ** Compiling $@\n"

libxdrfile.a:
	@cp $(SRCDIR)/libxdrfile.a .

xdr.o: $(SRCDIR)/xdr.f90
	@$(FC) -c $(SRCDIR)/xdr.f90 -cpp -g -fbacktrace 
	@printf "   ** Compiling $@\n"

mod_loops.o: $(SRCDIR)/mod_loops.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90 \
	$(SRCDIR)/mod_pot_ljc.f90 \
	$(SRCDIR)/mod_pot_bhc.f90 \
	$(SRCDIR)/mod_pot_ljc_pair.f90 \
	$(SRCDIR)/xdr.f90 
	@$(FC) -c $(SRCDIR)/mod_loops.f90 $(FFLAGS) 
	@printf "   ** Compiling $@\n"

#resume

mod_write_vmd.o:
	@$(FC) -c $(SRCDIR)/mod_write_vmd.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_search_structures.o: $(SRCDIR)/mod_search_structures.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90 \
	$(SRCDIR)/mod_loops.f90 
	@$(FC) -c $(SRCDIR)/mod_search_structures.f90 $(FFLAGS) 
	@printf "   ** Compiling $@\n"

mod_resume.o: $(SRCDIR)/mod_resume.f90 \
	$(SRCDIR)/mod_cmd_line.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_grids.f90 \
	$(SRCDIR)/mod_search_structures.f90
	@$(FC) -c $(SRCDIR)/mod_resume.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_deallocate_all.o: $(SRCDIR)/mod_deallocate_all.f90 \
	$(SRCDIR)/mod_cmd_line.f90 \
	$(SRCDIR)/mod_read_molecules.f90 \
	$(SRCDIR)/mod_grids.f90 \
	$(SRCDIR)/mod_input_read.f90 \
	$(SRCDIR)/mod_pot_ljc.f90 \
	$(SRCDIR)/mod_pot_bhc.f90 \
	$(SRCDIR)/mod_pot_ljc_pair.f90 \
	$(SRCDIR)/mod_search_structures.f90
	@$(FC) -c $(SRCDIR)/mod_deallocate_all.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

clean: 
	@printf "\n"
	@printf "   ** Deleting *.o files\n"
	@rm -rf $(OBJS) 
	@printf "\n"
	@printf "   ** Deleting *.mod files\n"
	@rm -rf $(MODS) 
	@printf "\n"
	@printf "   ** Deleting $(PROGRAM) file\n"
	@rm -rf $(PROGRAM)
	@printf "\n"

############################################
## install: $(PROGRAM)										##
##	install -t /usr/local/bin $(PROGRAM)  ##
############################################
