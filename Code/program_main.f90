!##################################################!
!### Main Molecular Dynamics simulation program ###!
!##################################################!

program main
implicit none
integer :: M              !Number of particles per dimension in FCC lattice
integer :: N              !Total number of particles
integer :: nhis,ngr       !RDF subroutine counters
integer :: i,counter
real*8  :: upot           !Potential Energy (reduced units)
real*8  :: time           !Time (reduced units)
real*8  :: density        !Density of the lattice, reduced units
real*8  :: cut            !LJ forces calculation distance cut-off
real*8  :: L              !Simulation box width
real*8  :: dt             !Time step for integration (reduced units)
real*8  :: Temp           !Simulation temperature (reduced units)
real*8  :: press          !Term related to pressure calculated in the forces subroutine
real*8  :: kin            !Kinetic energy (reduced units)
real*8  :: delg           !Bin size in RDF calculation
real*8  :: mass           !Atomic mass (g/mol)
real*8  :: eps            !LJ epsilon parameter (J/mol)
real*8  :: sigma          !LJ sigma parameter (J/mol)
real*8  :: Maxtime        !Maximum simulation time
real*8,allocatable :: r(:,:)      !Matrix with the xyz coordinates of all particles at a given step
                                  !(reduced units)
real*8,allocatable :: vel(:,:)    !Matrix with the velocities of all particles at a given step
                                  !(Reduced units)
real*8,allocatable :: force(:,:)  !Matrix with forces at a given step (Reduced units)
real*8,allocatable :: g(:)        !Radial Distribution function


!### Read input parameters from input.dat ###!

call input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime) 


!### With the M given, calculate the total number of particles ###!

N=M**(3)*4

!### Allocate matrixes and vectors ###!

allocate(r(N,3),vel(N,3),force(N,3),g(nhis))

!### System initialization ###!

!Generate FCC lattice
call coordenadas(N,M,r,density,L)

print *, L

!Adjust lattice with PBC
call boundary_conditions(r,N,L)

!Print initial configuration
call trajectory(r,N,time,counter)

!Assing initial velocities consistent with temperature
call in_velocity(vel,N,Temp)

!Initialize RDF
call rdf(r,N,L,0,nhis,density,delg,ngr,g)

!Cut-off calculation according to box width
cut=L*0.5d0

!Time initialization
time=0.d0

!Calculate initial forces,energies,pressure
call forces_LJ(L,N,r,cut,force,press,upot)

!Print first values of the magnitudes
call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)

!### Main Molecular Dynamics loop ###!

do while (time.lt.Maxtime)
    time=time+dt
    counter=counter+1
    !New positions and velocities
    call verlet_velocity(N,cut,press,r,vel,force,dt,upot,L)

    if (time.gt.0.3*Maxtime)then !After it equilibrates

        !Kinetical energy calculation
        call kinetic_en(vel,N,kin)
        !Instant temperature calculation
        call temperatura(kin,N,Temp)
        !Print positions
        call trajectory(r,N,time,counter)
        !Pirnt magnitudes
        call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)
        !Update RDF
        call rdf(r,N,L,1,nhis,density,delg,ngr,g)
        !Print total momentum
        call momentum(time,vel,N)

    end if
end do

!Final RDF calculation
call rdf(r,N,L,2,nhis,density,delg,ngr,g)

end program main
    
