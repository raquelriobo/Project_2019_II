! ** Kinetic energy subroutine ** !
subroutine kinetic_en(vel,N,kin)
implicit none
real*8 :: kin,vel(N,3)
integer :: i,N
kin=0.0
do i=1,N
  kin=kin+0.5d0*(vel(i,1)**2d0+vel(i,2)**2d0+vel(i,3)**2d0)
end do
endsubroutine kinetic_en
