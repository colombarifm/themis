# makefile for themis program

# program name
PROGRAM = com


# version control 
# Adding revision informatio to themis

#TODO:
#REV = $(shell git rev-parse --short HEAD)
# an error occurs

FC = gfortran

#FFLAGS = -g -Wall -O3 -ffpe-trap=invalid,zero,overflow -fcheck=bounds -Wno-compare-real -Wno-conversion -fbacktrace -fcheck=all -Wextra 
FFLAGS = -g -Wall -O3 -fcheck=bounds -Wno-compare-reals -Wno-conversion -fbacktrace -fcheck=all -Wextra 

# source directory
SRCDIR = src

# source files and objects
SRC = $(wildcard $(SRCDIR)/com.f90)

# objects created
OBJS = mod_error_handling.o \
	mod_constants.o \
  mod_info.o \
	mod_cmd_line.o \
	mod_inquire.o \
	mod_read_molecule.o \
	mod_calc_com.o \
	mod_deallocate_all.o

# modules created
MODS = mod_error_handling.mod \
	mod_constants.mod \
  mod_info.mod \
	mod_constants.mod \
	mod_cmd_line.mod \
	mod_inquire.mod \
	mod_read_molecule.mod \
	mod_calc_com.mod \
	mod_deallocate_all.mod

all: $(OBJS)
	@printf "\n"
	@printf "\n"
	@printf "  GENERATING EXECUTABLE...\n"
	@printf "\n"
	@$(FC) $(SRC) $(FFLAGS) $(OBJS) -o $(PROGRAM) 
	@printf "   ** Compiling com\n"
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

mod_read_molecule.o: $(SRCDIR)/mod_read_molecule.f90 \
	$(SRCDIR)/mod_cmd_line.f90 
	@$(FC) -c $(SRCDIR)/mod_read_molecule.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_calc_com.o: $(SRCDIR)/mod_calc_com.f90 \
	$(SRCDIR)/mod_constants.f90 \
	$(SRCDIR)/mod_cmd_line.f90 \
	$(SRCDIR)/mod_read_molecule.f90
	@$(FC) -c $(SRCDIR)/mod_calc_com.f90 $(FFLAGS)
	@printf "   ** Compiling $@\n"

mod_deallocate_all.o: $(SRCDIR)/mod_deallocate_all.f90 \
	$(SRCDIR)/mod_read_molecule.f90
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
