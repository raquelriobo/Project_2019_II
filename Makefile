#Makefile


#Compilador
F90=gfortran

#Programa principal
TARGET=program_main

#Generación de plots
results.png: Results.txt
	gnuplot Scripts_GNUPlot/plot_Energy_Raquel.gnu
	gnuplot Scripts_GNUPLot/plot_Moment_Raquel.gnu

#Ejecución del programa
Results.txt : $(TARGET).x 
	./$(TARGET).x < input.txt

#Compilación del programa completo
$(TARGET).x : $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o fuerzas_Raquel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o

	$(F90) -o $(TARGET).x $(TARGET).o Ekinetic_Raquel.o boundary.o verlet_vel.o fuerzas_Raquel.o in_velocities_Raquel.o coordenadas_Raquel.o units_print.o temperatura.o trajectory.o radial.o moment_Raquel.o forces_RaquelNEW.o input.o


#Todos los archivos con extensión .f90 se compilan en archivos .o
%.o : %.f90
	$(F90) -c $<


##Comentarios durante la ejecución del Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile



##backup : Copia de seguridad comprimida
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90

##clean : Rule to clean executable objects 
.PHONY : clean
clean:
	rm *.o
