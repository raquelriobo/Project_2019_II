! ** Kinetic energy function ** !
subroutine kinetic_en(vel,N)
!real*8 function kinetic_en(vel,N)
implicit none
real*8 :: kin,vel(N,3)
integer :: i,N
kinetic_en=0.0
do i=1,N
  kinetic_en=kinetic_en+0.5d0*(vel(i,1)**2d0+vel(i,2)**2d0+vel(i,3)**2d0)
end do
!end function kinetic_en
endsubroutine kinetic_en
