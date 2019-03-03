! ### FORCE ### !
subroutine forces_LJ(L,N,r,cut,force,press,upot)
implicit none
integer :: i,j
integer :: N                    ! Number of particles
real(8) :: L                    ! Cell long
real(8) :: r(N,3)               ! Coordinates matrix
real(8) :: cut                  ! Cut-off
real(8) :: dx,dy,dz             ! Distance components
real(8) :: ff, d2, d2i, d6i
real(8) :: force(N,3)           ! Forces matrix
real(8) :: upot                 ! Potential energy
real(8) :: press                ! Pressure

force=0d0
upot=0d0
press=0d0

do i=1,N-1
    do j=1+i,N                  ! Loop over all pairs
        dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
        dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
        dy=r(i,2)-r(j,2)
        dy=dy-L*nint(dy/L)
        dz=r(i,3)-r(j,3)
        dz=dz-L*nint(dz/L)

        d2=dx**2d0+dy**2d0+dz**2d0

        if(d2.lt.cut**2)then    ! Test cutoff
        d2i=1/d2
        d6i=d2i**3

            ! ### Compute Forces ### !
            ff=(48d0*d2i*d6i*(d6i-0.5))
            force(i,1)=force(i,1)+ff*dx
            force(i,2)=force(i,2)+ff*dy
            force(i,3)=force(i,3)+ff*dz
            force(j,1)=force(j,1)-ff*dx
            force(j,2)=force(j,2)-ff*dy
            force(j,3)=force(j,3)-ff*dz

            ! ### Compute Potential Energy ### !
            upot=upot+4d0*d6i*(d6i-1)

            ! ### Compute Pressure ### !
            press=press+ff*dx*dx+ff*dy*dy+ff*dz*dz
        end if
    end do
end do
end subroutine forces_LJ
