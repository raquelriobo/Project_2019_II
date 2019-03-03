subroutine coordenadas(N,M,coord,density,L)
implicit none
integer :: i,j,k,cont
integer :: N				! Number of particles
integer :: M				!
real(8) :: density, L, a
real(8) :: r01(3),r02(3),r03(3),r04(3)	! Initial coordinates
real(8) :: vec(3)
real(8) :: coord(N,3)			! Total coordinates

! ### Data ### !
L=(N/density)**(1d0/3d0)		! Cell long
a=L/M 					! Distance between particles in the cell

! ### Face Centered Cubic Coordinates ### !
r01(1)=0.0;r01(2)=0.0;r01(3)=0.0
r02(1)=a/2;r02(2)=a/2;r02(3)=0.0
r03(1)=0.0;r03(2)=a/2;r03(3)=a/2
r04(1)=a/2;r04(2)=0.0;r04(3)=a/2

! ### Add Vectors ### !
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

end subroutine coordenadas
