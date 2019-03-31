##Makefile for a Molecular Dynamics simulation program

#**Variables**#

#Compiler (Serial code)
F90=gfortran

#Compiler (Parallel code)
MPI=mpif90

#Number of workers
PROC=2

#Main program
TARGET=program_main


#Energy, pressure, temperature and total momentum plot generation
energies.eps : Results.txt
	@echo "Generating plots with the results..."
	@gnuplot Scripts_GNUPlot/plot_Energy_Raquel.gnu
	@gnuplot Scripts_GNUPlot/plot_momentum.gnu
	@gnuplot Scripts_GNUPlot/plot_rdf.gnu
	@echo Copying to Results folder...
	@mkdir -p Results
	@cp *.txt *.eps *.xyz Results
	@echo "Done!"

#Main program execution
Results.txt : $(TARGET).x
	@echo "Executing the program with the input values and $(PROC) workers." 
	mpirun -np $(PROC) $(TARGET).x
	@echo Generating Results.txt file ...
	./Scripts_GNUPlot/fix_results.sh

#Compilation of the main program
$(TARGET).x : $(TARGET).o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o input.o kin_mpi_new.o mpi_initial_velocity.o boundary_mpi.o verlet_mpi.o verlet_list_mpi.o verlet_list.o forces_vlist_RaquelNEW.o checkverletlist.o radial_mpi.o mpi_moment.o
	@echo "Compiling program_main.x ..."
	$(MPI) -g -o $(TARGET).x $(TARGET).o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o input.o kin_mpi_new.o mpi_initial_velocity.o boundary_mpi.o verlet_mpi.o verlet_list_mpi.o verlet_list.o forces_vlist_RaquelNEW.o checkverletlist.o radial_mpi.o mpi_moment.o


#All files with extension .f90 are compiled to objects .o
%.o : Code_MPI/%.f90
	@echo "Compiling the necessary subroutines ..."
	$(MPI) -g -c $<


##statistics : binning of the time series for different magnitudes
.PHONY : statistics
statistics :
	@$(F90) -o binning.x Code_MPI/binning2.f90
	@echo Performing statistical analysis ...
	@./binning.x
	@gnuplot Scripts_GNUPlot/plot_binning.gnu
	@echo Moving results to Results_binning
	@mkdir -p Results_binning
	@cp *binning.txt *mit.txt *bin.eps Results_binning
	@echo Done!

##help: instructions about the use of this Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile

##backup : make a compressed copy of the base code
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" Code_MPI/*.f90

##clean_all : rule to clean executable objects, results and images
.PHONY : clean_all
clean_all:
	@echo Removing compiled objects and results
	rm *.o *.txt *.eps *.xyz *.x

##clean_exe : rule to clean only executable objects
.PHONY : clean_exe
clean_exe:
	@echo Removing compiled objects only
	rm *.o *.x

##variables : print variables of this Makefile
.PHONY : variables
variables :
	@echo List of variables:
	@echo Compiler serial: F90=$(F90)
	@echo Compiler parallel: MPI=$(MPI)
	@echo Number of workers: PROC=$(PROC)
	@echo Main target: TARGET=$(TARGET)

##Nota sobre compilación: si se modifica el código fuente de una subrutina (no program_main), se puede usar la rule clean_exe y después make, o se puede usar "make nombresubrutina.o" para compilar de nuevo esa subrutina, y después make program_main.xpara recompilar el ejecutable final.
