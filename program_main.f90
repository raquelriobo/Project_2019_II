!Programa Principal


program main
implicit none
integer :: M,N,Maxstep
real*8 :: density,cut,L,dt,Temp
real*8,allocatable :: r(:,:),vel(:,:),force(:,:)

!Parameters
M=10
density=0.9
dt=0.0003
mass=4 
eps=91.04  !J/mol
sigma=2.963d-10 !m


N=M**3.d0


allocate(r(N,3),vel(N,3),force(N,3))

call coordenadas_raquel(M,N,r,density,L)

!possible pbc

Temp=3d2 !Hay q pasarlos a unidades reducidas

call in_velocities(vel,N,Temp)

call verlet_velocity(N,r,vel,mass,eps,sigma)

end program main
