!##################################################!

program main
implicit none
integer :: M              !Number of particles per dimension in FCC lattice
integer :: N              !Total number of particles
integer :: nhis,ngr       !RDF subroutine counters
integer :: i,counter,cont,contF
integer, allocatable :: nlist(:),list(:,:)
real*8  :: upot           !Potential Energy (reduced units)
real*8  :: time           !Time (reduced units)
real*8  :: density        !Density of the lattice, reduced units
real*8  :: cut            !LJ forces calculation distance cut-off
real*8  :: cut2            !LJ forces calculation distance cut-off
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
real*8  :: resta
real*8  :: start, finish
real(8) :: suma,media
real(8) :: timeforces(150000)
real*8,allocatable :: r(:,:)      !Matrix with the xyz coordinates of all particles at a given step
real*8,allocatable :: rold(:,:)      !Matrix with the xyz coordinates of all particles at a given step
                                  !(reduced units)
real(8) :: dx,dy,dz
real*8,allocatable :: vel(:,:)    !Matrix with the velocities of all particles at a given step
                                  !(Reduced units)
real*8,allocatable :: force(:,:)  !Matrix with forces at a given step (Reduced units)
real*8,allocatable :: g(:)        !Radial Distribution function
contF=1

call cpu_time(start)

!### Read input parameters from input.dat ###!
call input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime) 


!### With the M given, calculate the total number of particles ###!
N=M**(3)*4

!### Allocate matrixes and vectors ###!
allocate(r(N,3),vel(N,3),force(N,3),g(nhis),rold(N,3))
allocate(nlist(N),list(N,N-1))

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
cut=L*0.4d0
!cut=L*0.3d0
cut2=L*0.5d0
!Time initialization
time=0.d0
counter=0
cont=0
!Calculate initial forces,energies,pressure
call new_vlist(i,L,N,r,cut2,nlist,list,rold)

call forces_vlist(L,N,r,cut,force,press,upot,nlist,list,contF,timeforces)

!### Main Molecular Dynamics loop ###!

do while (time.lt.Maxtime)
    time=time+dt
    counter=counter+1
    
!New positions and velocities
call verlet_velocity(N,cut,press,r,vel,force,dt,upot,L,nlist,list,contF,timeforces)

!Print positions
if (mod(counter,100).eq.0) then
    call trajectory(r,N,time-time*0.3,counter)
end if

if (mod(counter,500).eq.0) then
    call checkvlist(i,L,N,r,cut,cut2,nlist,list,rold,dx,dy,dz,counter,cont)    
endif

!Kinetical energy calculation
    call kinetic_en(vel,N,kin)
!Instant temperature calculation
    call temperatura(kin,N,Temp)
!Print total momentum
    call momentum(time-time*0.3,vel,N)
!Pirnt magnitudes
    call units_print(time-time*0.3,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)

    if (time.gt.0.4*Maxtime)then !After it equilibrates measure the properties

!Update RDF
        if (mod(counter,100).eq.0) then
            call rdf(r,N,L,1,nhis,density,delg,ngr,g)
        end if
    end if
end do

! ### Final RDF calculation
call rdf(r,N,L,2,nhis,density,delg,ngr,g)

! ### Calculo CPU time ### !
    suma=sum(timeforces(:contF))
    media=suma/contF
    print*, 'Suma Fuerzas', suma
    print*, 'Media Fuerzas', media

call cpu_time(finish)
print*, 'Total time=',finish-start, 'seconds'
end program main
    
