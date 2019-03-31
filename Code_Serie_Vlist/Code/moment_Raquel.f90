! ### Total Momentum ### !
subroutine momentum(time,vel,N)
implicit none
integer :: i
integer :: N        ! Number of particles
real(8) :: p(3)     ! Total Momentum vector
real(8) :: vel(N,3) ! Velocity matrix
real(8) :: time

open(70,file="Momentum.txt")

p=0d0

do i=1,N                ! Loop over particles
    p(1)=p(1)+vel(i,1)  ! px
    p(2)=p(2)+vel(i,2)  ! py
    p(3)=p(3)+vel(i,3)  ! pz
end do

write(70,*) time,p(1),p(2),p(3)

endsubroutine momentum

