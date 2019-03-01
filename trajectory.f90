
subroutine trajectory(r,N,time)
implicit none
integer :: N,i
real*8 :: r(N,3)
real*8 :: time

open(14,file="TRAJ.xyz")

write(14,*)N
write(14,*) "time=",time
do i=1,N
  write(14,*)"He",r(i,1),r(i,2),r(i,3)
end do

end subroutine trajectory
