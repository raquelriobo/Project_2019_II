! ***** Velocities inialization **** !
subroutine in_velocity(vel,N,Temp)
implicit none
integer :: N,i,j
real*8 :: Temp,kin,vel(N,3),suma(N),kinetic_en

! ** Random velocities ** !
do i=1,N
  vel(i,:)=(2d0*rand()-1d0)/2d0
end do

! ** Kinetic energy ** !
kin=kinetic_en(vel,N)

! ** Reescaling Kinetic energy ** !
vel=vel*sqrt(N*3d0*Temp/(2d0*kin))

! ** Total velocity = 0 ** !
suma=0d0
do i=1,3
  do j=1,N
    suma(i)=suma(i)+vel(j,i)
  end do
end do
suma=suma/N
do i=1,N
  do j=1,3
    vel(i,j)=vel(i,j)-suma(j)
  end do
end do
end subroutine in_velocity

! ** Kinetic energy function ** !
real*8 function kinetic_en(vel,N)
implicit none
real*8 :: kin,vel(N,3)
integer :: i,N
kinetic_en=0.0
do i=1,N
  kinetic_en=kinetic_en+0.5d0*(vel(i,1)**2d0+vel(i,2)**2d0+vel(i,3)**2d0)
end do
end function kinetic_en
