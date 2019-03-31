!###########################!
!### Trajectory printing ###!
!###########################!

subroutine trajectory(r,N,time,counter)
implicit none
integer :: N         !Number of particles
integer :: i,counter
real*8 :: r(N,3)     !Coordinates matrix  
real*8 :: time       !Time (reduced)

open(14,file="TRAJ.xyz")

if (mod(counter,100).eq.0) then !Print every 100
    write(14,*)N
    write(14,*) "time=",time
    do i=1,N
        write(14,*)"He",r(i,1),r(i,2),r(i,3)
    end do
end if
end subroutine trajectory
