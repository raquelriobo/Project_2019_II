!Programa Principal


program main
implicit none
integer :: N,Maxstep
real*8 :: density,cut,L,dt,Temp
real*8,allocatable :: r(:,:),vel(:,:),force(:,:)

!Parameters
N=108
density=0.9
dt=0.0003

allocate(r(N,3),vel(N,3),force(N,3))


