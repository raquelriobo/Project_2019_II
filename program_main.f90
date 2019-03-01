!Programa Principal


program main
implicit none
integer :: M,N,i
real*8 :: upot,time,density,cut,L,dt,Temp,press,kin
real*8 :: mass,eps,sigma,Maxtime
real*8,allocatable :: r(:,:),vel(:,:),force(:,:)
real*8 :: kb=1.3806d-23 !J/K

!Parameters
M=3
density=0.9d0 !Reduced
dt=0.0003d0 !Reduced
mass=4.0d0 !g/mol 
!eps=91.04d0  !J/mol
eps=0.998d3
!sigma=2.963d-10 !m
sigma=3.4d-10
!Temp=(kb*300.0d0)/(eps*6.022d-23) !Reduced units
Temp=10
Maxtime=1
N=M**(3)*4 !Number of particles

allocate(r(N,3),vel(N,3),force(N,3))


!System initialization

call coordenadas(N,M,r,density,L)

!Cut-off calculation

cut=L*0.3d0
time=0.d0

call boundary_conditions(r,N,L)

call trajectory(r,N,time)

call in_velocity(vel,N,Temp)

call forces_LJ_Press(L,N,r,cut,force,press,upot)

do while (time.lt.Maxtime)
  call verlet_velocity(N,cut,press,r,vel,force,dt,upot)
  call boundary_conditions(r,N,L)
  call kinetic_en(vel,N,kin)
  call temperatura(kin,N,Temp)
  call trajectory(r,N,time)
  call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)
  time=time+dt
end do


end program main
