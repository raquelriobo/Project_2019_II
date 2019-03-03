
subroutine trajectory(r,N,time,counter)
implicit none
integer :: N,i,counter
real*8 :: r(N,3)
real*8 :: time
open(14,file="TRAJ.xyz")

if (mod(counter,100).eq.0) then
    write(14,*)N
    write(14,*) "time=",time
    do i=1,N
        write(14,*)"He",r(i,1),r(i,2),r(i,3)
    end do
end if
end subroutine trajectory
