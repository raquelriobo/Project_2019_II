# Project_2019_II

Molecular Dynamics simulation program

Authors:

- Raquel Agraso:
	Initial coordinates (coordenadas_Raquel.f90)
	Initial velocities (in_velocities_Raquel.f90)
	Forces/Potential Energy calculation (forces_RaquelNEW.f90)
	Total Momentum calculation (moment_Raquel.f90)
	Kinetic energy calculation (Ekinetic_Raquel.f90)
	GNUPLot scripts
- Beatriz Piniello:
	Main program (program_main.f90)
	Makefile
	Printing trajectory for visualization (trajectory.f90)
	Radial Distribution Function (radial.f90)
	GNUPlot scripts
- Encarna Vallès:
	Periodic Boundary Conditions (boundary.f90)
	Reading input parameters (input.f90)
	Instant temperature calculation (temperatura.f90)
	Velocity verlet integrator (verlet_vel.f90)
- Alba Villar:
	Printing subroutine (units_print.f90)
	Binning subroutine (binning2.f90)
	GNUPlot scripts

General repository information:

Folders:

- Code: all necessary .f90 files
- Inputs: different input files for calculation
- Scripts_GNUPlot: scripts to obtain figures

Execution instructions:

- To obtain results with the inputs found in input.dat, just type 'make' in the main directory
- To obtain the binning plots, 'make statistics'
- For more options in Makefile, 'make help'

Files:

1. coordenadas_Raquel.f90   : Initiallice positions with a FCC structure according with the pressure

2. in_velocities_Raquel.f90 : Intitiallice velocities

3. forces_RaquelNEW.f90     : Compute the forces that act on each particle

4. moment_Raquel.f90        : Compute the different components of the total momentum

5. Ekinetic_Raquel.f90      : Compute the Kinetic energy

6. trajectory.f90           : 

7. radial.f90               :

8. boundary.f90             :

9. temperatura.f90          :

10. verlet_vel.f90          :

11. units_print.f90         :

12. binning2.f90            : Compute the

13. program_main.f90        :

14. Makefile                :

15. input.f90               :

6. GNUPLot scripts          :

