! ### Kinetic energy ### !
subroutine kinetic_en(vel,N,kin)
implicit none
integer :: i
integer :: N            ! Number of particles
real(8) :: kin          ! Kinetic energy
real(8) :: vel(N,3)     ! Velocity matrix

kin=0.0

do i=1,N
    kin=kin+0.5d0*(vel(i,1)**2d0+vel(i,2)**2d0+vel(i,3)**2d0)
end do
endsubroutine kinetic_en
