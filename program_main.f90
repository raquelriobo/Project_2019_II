!Programa Principal


program main
implicit none
integer :: M,N,Maxstep,i
real*8 :: upot,time,density,cut,L,dt,Temp,press,kin
real*8 :: mass,eps,sigma
real*8,allocatable :: r(:,:),vel(:,:),force(:,:)
real*8 :: kb=1.3806d-23 !J/K

!Parameters
M=10
density=0.9d0 !Reduced
dt=0.0003d0 !Reduced
mass=4.0d0 !g/mol 
eps=91.04d0  !J/mol
sigma=2.963d-10 !m
Temp=(kb*300.0d0)/eps !Reduced units
MaxStep=1000
N=M**3 !Number of particles


allocate(r(N,3),vel(N,3),force(N,3))


!System initialization

call coordenadas(N,M,r,density,L)

!Cut-off calculation

cut=L*0.3d0

!(Possible pbc call)

call in_velocity(vel,N,Temp)

do i=1,MaxStep
!  call verlet_velocity(N,cut,press,r,vel,dt,upot)
!  call boundary_conditions(r,N,L)
!  call kinetic_en(vel,N,kin)
!  call temperatura(kin,N,Temp)
!  call units_print(time,upot,kin,press,L,dt,sigma,eps,density)
end do


end program main
