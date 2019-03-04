##Makefile for a Molecular Dynamics simulation program


#Compiler
F90=gfortran

#Main program
TARGET=program_main


#Energy, pressure, temperature and total momentum plot generation
energy.eps : Results.txt
	@echo "Generating plots with the results..."
	gnuplot Scripts_GNUPlot/plot_Energy_Raquel.gnu
	@echo "Done!"
total_momentum.eps : Momentum.txt
	gnuplot plot_momentum.gnu
rdf.eps : radial.txt
	gnuplot Scripts_GNUPlot/plot_rdf.gnu

#Main program execution
Results.txt : $(TARGET).x Inputs/input.dat
	@echo "Executing the program with the input values ..." 
	./$(TARGET).x < Inputs/input.dat

#Compilation of the main program
$(TARGET).x : $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o
	@echo "Compiling program_main.x ..."
	$(F90) -o $(TARGET).x $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o


#All files with extension .f90 are compiled to objects .o
%.o : Code/%.f90
	@echo "Compiling the necessary subroutines ..."
	$(F90) -c $<


##statistics : binning of the time series for different magnitudes
.PHONY : statistics
statistics :
	$(F90) -o binning.x binning2.f90
	./binning.x
	gnuplot Scripts_GNUPlot/plot_binning.gnu

##help: instructions about the use of this Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile


##backup : make a compressed copy of the base code
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90

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
