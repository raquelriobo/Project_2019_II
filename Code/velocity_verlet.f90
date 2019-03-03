subroutine velocity_verlet(N,cut,press,r,vel,force,dt,upot,L)
implicit none
real*8 :: r(N,3),vel(N,3),force(N,3)
real*8 :: dt
integer :: i,j,N
real*8 :: L,cut,upot,press

do j=1,N !Loop for the particles
      !New positions for each dimension
      r(j,1)=r(j,1)+vel(j,1)*dt+0.5d0*force(j,1)*dt**2d0
      r(j,2)=r(j,2)+vel(j,2)*dt+0.5d0*force(j,2)*dt**2d0
      r(j,3)=r(j,3)+vel(j,3)*dt+0.5d0*force(j,3)*dt**2d0

      !Velocity calculation (first part, old forces)
      vel(j,1)=vel(j,1)+0.5d0*force(j,1)*dt
      vel(j,2)=vel(j,2)+0.5d0*force(j,2)*dt
      vel(j,3)=vel(j,3)+0.5d0*force(j,3)*dt
    end do
    !Force and potential energy calculation in new time
    call forces_LJ_Press(L,N,r,cut,force,press,upot)
    do j=1,N
      !Velocity calculation (second part, new forces)
      vel(j,1)=vel(j,1)+0.5d0*force(j,1)*dt
      vel(j,2)=vel(j,2)+0.5d0*force(j,2)*dt
      vel(j,3)=vel(j,3)+0.5d0*force(j,3)*dt
    end do

    call boundary_conditions(r,N,L)
     
!    r=r-nint(r/L)*L
end subroutine velocity_verlet
