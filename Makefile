#Makefile


#Compilador
F90=gfortran

#Programa principal
TARGET=program_main

#Energy, pressure, temperature and total momentum plot generation
results.png: Results.txt
	@echo "Generating plots with the results..."
	gnuplot Scripts_GNUPlot/plot_Energy_Raquel.gnu
	gnuplot Scripts_GNUPLot/plot_Moment_Raquel.gnu
	@echo "Done."

#Main program execution
Results.txt : $(TARGET).x
	@echo "Executing the program with the input values ..." 
	./$(TARGET).x < input.dat

#Compilation of the main program
$(TARGET).x : $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o fuerzas_Raquel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o
	@echo "Compiling program_main.x ..."
	$(F90) -o $(TARGET).x $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o fuerzas_Raquel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o


#All files with extension .f90 are compiled to objects .o
%.o : %.f90
	@echo "Compiling the necessary subroutines ..."
	$(F90) -c $<


##help: instructions about the use of this Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile



##backup : make a compressed copy of the base code
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90

##clean : rule to clean executable objects and txt files
.PHONY : clean
clean:
	@echo Removing comiled objects and results
	rm *.o *.txt *.png *.xyz *.x
