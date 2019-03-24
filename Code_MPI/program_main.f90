!##################################################!
!### Main Molecular Dynamics simulation program ###!
!##################################################!

program main
implicit none
include 'mpif.h'
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

integer :: ierr
integer :: rank                     !Task name
integer :: size                     !Total number of processors
integer :: root                     !Processor 0
integer :: part1                    !Main division of particles/processor
integer :: part2                    !Uneven division of particles/processor
integer :: nini,nfin,nini_first,nfin_first !Indexes for matrix division
integer :: sizes(2),subsizes(2),starts(2)
integer :: blocktype,resizedtype
integer :: intsize,N_aux
integer (kind=MPI_Address_kind) :: start, extent


root=0
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)


!### Read input parameters from input.dat ###!

!call input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime) 

M=3
dt=0.0003
density=0.9
Temp=10
sigma=3.4d-10
eps=0.9d3
nhis=1000
Maxtime=1
mass=40

!### With the M given, calculate the total number of particles ###!

N=M**(3)*4

!### Division of particles among processors ###!

if (mod(N,size).eq.0)then
  part1=N/(size)
  part2=part1
else
  N_aux = N - mod(N,size)
  part1=int(N_aux/size)
  part2= part1 + mod(N,size)
end if

!### Indexes for matrix division ###!

nini_first=1
nfin_first=part2

nini=nfin_first+(rank-1)*part1
nfin=N


!### Subarray type creation ###!

sizes=  [part1*size, 3 ] !sizes array total
subsizes = [ part1, 3 ] !subsizes subarrays
starts   = [ 0, 0 ]


call MPI_Type_create_subarray( 2, sizes, subsizes, starts,     &
                                 MPI_ORDER_FORTRAN, MPI_REAL8, &
                                 blocktype, ierr)
start = 0
call MPI_Type_size(MPI_REAL8, intsize, ierr)
extent = intsize * part1

call MPI_Type_create_resized(blocktype, start, extent, resizedtype, ierr)
call MPI_Type_commit(resizedtype, ierr)


!### Allocate matrixes and vectors ###!

allocate(r(N,3),vel(N,3),force(N,3),g(nhis))

!### System initialization ###!

!Generate FCC lattice
call coordenadas(N,M,r,density,L)


!Adjust lattice with PBC
call boundary_conditions(r,N,L)


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

!### Main Molecular Dynamics loop ###!

do while (time.lt.Maxtime)
    time=time+dt
    counter=counter+1
    !New positions and velocities
    call verlet_mpi(N,cut,press,r,vel,force,dt,upot,L,root,part1,rank)
!    call verlet_velocity(N,cut,press,r,vel,force,dt,upot,L)
   
    if (time.gt.0.4*Maxtime)then !After it equilibrates measure the properties

        !Kinetical energy calculation
        call ekin_mpi(vel,N,kin,part1,root,rank)
        !Instant temperature calculation
        call temperatura(kin,N,Temp)
        !Print positions
        call trajectory(r,N,time,counter)
        !Print magnitudes
        call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)
        !Update RDF
        call rdf(r,N,L,1,nhis,density,delg,ngr,g)
        !Print total momentum
        call momentum(time,vel,N)

    end if
end do

!Final RDF calculation
call rdf(r,N,L,2,nhis,density,delg,ngr,g)



call MPI_FINALIZE(ierr)

end program main
    
