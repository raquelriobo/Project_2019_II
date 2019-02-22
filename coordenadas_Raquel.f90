subroutine coordenadas(N,coord,density,L)
implicit none
integer :: N,M,i,j,k,cont
real*8 :: a,L,r01(3),r02(3),r03(3),r04(3),vec(3),coord(N,3),density
L=(N/density)**(1d0/3d0)  !Cell long
!M=(N/4d0)**(1d0/3d0) !Number of times fcc is replicated in each dimension
N=(Md3)*4
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
