program Dynamics_Final2
implicit none
integer :: N,Maxstep,hh,i
real*8 :: density,cut,L,sigma,eps,d,dt,upot,kin,Temp
real*8,allocatable :: coord(:,:),vel(:,:),force(:,:),press
character(20) :: file_name

! ********** DATA ********** !
N=108 !Nº particles
dt=0.0003 !Time step
cut=2.0 !Truncation

allocate(coord(N,3),vel(N,3),force(N,3))
density=0
DO hh=0,3
density=density+0.2 ! **Density (reduced units)
write(6,*) density

! ******** COORDENADAS FACE CENTERED CUBIC ******** !
call coordenadas(N,coord,density,L)

! ** Initial temperature and steps to 'melt' fcc lattice
Temp=1000d0
MaxStep=2000d0

! ******** Velocity inialization ******** !
call in_velocity(vel,N,Temp)

file_name="v_verlet_1_"

END DO
end program Dynamics_Final2

subroutine coordenadas(N,coord,density,L)
implicit none
integer :: N,M,i,j,k,cont
real*8 :: a,L,r01(3),r02(3),r03(3),r04(3),vec(3),coord(N,3),density
L=(N/density)**(1d0/3d0)  !Cell long
M=(N/4d0)**(1d0/3d0) !Number of times fcc is replicated in each dimension
a=L/M !Distance between particles in the cell

! ********** COORDENADAS FACE CENTERED CUBIC********** !
r01(1)=0.0;r01(2)=0.0;r01(3)=0.0
r02(1)=a/2;r02(2)=a/2;r02(3)=0.0
r03(1)=0.0;r03(2)=a/2;r03(3)=a/2
r04(1)=a/2;r04(2)=0.0;r04(3)=a/2

! ********** Suma vectores ********** !
cont=1
do i=1,M
  do j=1,M
    do k=1,M
      vec=(/i,j,k/)
      coord(cont,:)=r01(:)+a*vec(:)
      coord(cont+1,:)=r02(:)+a*vec(:)
      coord(cont+2,:)=r03(:)+a*vec(:)
      coord(cont+3,:)=r04(:)+a*vec(:)
      cont=cont+4
    end do
  end do
end do

! ** Periodic Boundary Conditions
coord=coord-nint(coord/L)*L !If particle is further than L/2

end subroutine coordenadas

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

