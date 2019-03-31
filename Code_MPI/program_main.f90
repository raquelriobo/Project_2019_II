!#################################################!
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
real*8  :: cut2
real*8  :: resta
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
real*8,allocatable :: rold(:,:)
real*8,allocatable :: vel(:,:)    !Matrix with the velocities of all particles at a given step
                                  !(Reduced units)
real*8,allocatable :: force(:,:)  !Matrix with forces at a given step (Reduced units)
real*8,allocatable :: g(:)        !Radial Distribution function
integer, allocatable :: nlist(:),list(:,:) !Vector of Verlet list longitudes, and Verlet lists matrix 

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
logical :: check
real*8 :: forces_start(150000),forces_finish(150000),total_forces_time(150000)
real*8 :: start_time,finish_time,vlist_start,vlist_finish


root=0
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

!Measure time for the entire run

call cpu_time(start_time)

!### Read input parameters from input.dat ###!

call input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime,rank) 

!### With the M given, calculate the total number of particles ###!

N=M**(3)*4

!### Division of particles among processors ###!

if (size.gt.N) then
    write(6,*) "The number of workers cannot be higher than the number of particles"
    call MPI_ABORT(ierr)
    stop
end if


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

sizes=  [part1*size,3] !sizes array total
subsizes = [part1,3] !subsizes subarrays
starts   = [0,0]


call MPI_Type_create_subarray( 2, sizes, subsizes, starts,     &
                                 MPI_ORDER_FORTRAN, MPI_REAL8, &
                                 blocktype, ierr)
start = 0
call MPI_Type_size(MPI_REAL8, intsize, ierr)
extent = intsize * part1

call MPI_Type_create_resized(blocktype, start, extent, resizedtype, ierr)
call MPI_Type_commit(resizedtype, ierr)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!### Allocate matrixes and vectors ###!

allocate(r(N,3),vel(N,3),force(N,3),g(nhis),rold(N,3))
allocate(nlist(N),list(N,N-1))


!### System initialization ###!


!Generate FCC lattice
call coordenadas(N,M,r,density,L)


!Adjust lattice with PBC
call boundary_conditions_mpi(r,N,L,part1,part2,root,&
rank,nini,resizedtype,size)


!Assing initial velocities consistent with temperature
call in_velocity_mpi(vel,N,Temp,part1,part2,nini,nfin,&
root,rank,resizedtype,size)

!Initialize RDF
call rdf_mpi(r,N,L,0,nhis,density,delg,ngr,g,rank,root,nini,part1,part2,nlist,list)

!Cut-off calculation according to box width
cut=L*0.4d0   !Lennard-Jones cut-off
cut2=L*0.5d0  !Verlet list cut-off

!Time initialization
time=0.d0
counter=0
check=.false.

call cpu_time(vlist_start)

!Calculate initial Verlet lists
call vlist(i,L,N,r,cut2,nlist,list,rold)

call cpu_time(vlist_finish)

!Calculate initial forces,energies,pressure
call forces_vlist(L,N,r,cut,force,press,upot,nlist,list,rank,root,part1,part2,nini,size&
,resizedtype)

!### Main Molecular Dynamics loop ###!
do while (time.lt.Maxtime)
    time=time+dt
    counter=counter+1

    !New positions and velocities
    call verlet_mpi(N,cut,press,r,vel,force,dt,upot,L,root,rank,part1,part2,&
    nini,size,resizedtype,nlist,list,counter,forces_start,forces_finish)
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)   

    !Print positions
    if (mod(counter,100).eq.0) then
        call trajectory(r,N,time,counter,rank)
    end if

    if ((mod(counter,500).eq.0).and.rank.eq.root) then
        call checkvlist(i,L,N,r,cut,cut2,nlist,list,rold,counter,check)
    endif

        call MPI_BCAST(check,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

    if (check.eqv..true.) then
        call vlist(i,L,N,r,cut2,nlist,list,rold)
        check=.false.
    end if

        !Kinetical energy calculation
        call kinetic_en(vel,N,kin,part1,part2,nini,nfin,nini_first,nfin_first,&
        root,rank)

        !Instant temperature calculation
        call temperatura(kin,N,Temp)

        !Print total momentum
        call moment_mpi(time, vel, N, part1, part2, nini, &
        root, rank)

        !Print magnitudes
        call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass,rank)


    if (time.gt.0.4*Maxtime)then !After it equilibrates measure RDF

        !Update RDF
        if (mod(counter,100).eq.0) then
            call rdf_mpi(r,N,L,1,nhis,density,delg,ngr,g,rank,root,nini,part1,part2,nlist,list)
        end if

    end if
end do

!Final RDF calculation

call rdf_mpi(r,N,L,2,nhis,density,delg,ngr,g,rank,root,nini,part1,part2,nlist,list)

!Measure the time of the entire run
call cpu_time(finish_time)

!Time results:

if (rank.eq.root) then
    !write(6,*) "Total time:",finish_time-start_time, "seconds."
    !write(6,*) finish_time-start_time
    !write(6,*) "Total Verlet list time:", vlist_finish-vlist_start,"seconds."
    !write(6,*) vlist_finish-vlist_start
    total_forces_time=forces_finish-forces_start
    !write(6,*) "Total Force time (all iterations):",sum(total_forces_time),"seconds."
    !write(6,*) sum(total_forces_time)
    !write(6,*) "Average force time per iteration:",sum(total_forces_time)/counter,"seconds." 
    !write(6,*) sum(total_forces_time)/counter
     write(6,*) finish_time-start_time,vlist_finish-vlist_start,sum(total_forces_time),sum(total_forces_time)/counter
end if


call MPI_FINALIZE(ierr)

end program main
    
