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
