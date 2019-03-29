program main
implicit none
integer,parameter :: N=108
real(8) :: vel(N,3)
real(8) :: kin
real(8) :: cut                  ! Cut-off
real(8) :: force(N,3)
real(8) :: L
real(8) :: upot                 ! Potential energy
real(8) :: press 
real(8) :: r(N,3) 
!integer :: i,ierr,rank,size,root,part1

call forces_vlist(L,N,r,cut,force,press,upot)

end program main
