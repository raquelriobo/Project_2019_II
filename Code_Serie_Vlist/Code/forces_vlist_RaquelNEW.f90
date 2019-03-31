! ### FORCE ### !
subroutine forces_vlist(L,N,r,cut,force,press,upot,nlist,list,cont,timeforces)
implicit none
integer :: i,j,jj,cont
integer :: N                    ! Number of particles
integer :: nlist(N),list(N,N-1)
real(8) :: L                    ! Cell long
real(8) :: r(N,3)               ! Coordinates matrix
real(8) :: cut                  ! Cut-off
real(8) :: dx,dy,dz             ! Distance components
real(8) :: ff, d2, d2i, d6i, d
real(8) :: force(N,3)           ! Forces matrix
real(8) :: upot                 ! Potential energy
real(8) :: press                ! Pressure
real(8) :: startF(150000), finishF(150000),timeforces(150000)
real(8) :: media
call cpu_time(startF(cont))

force=0d0
upot=0d0
press=0d0
do i=1,N
  do jj=1,nlist(i)
    j=list(i,jj)
    if (j.gt.i) then
            dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
            dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
            dy=r(i,2)-r(j,2)
            dy=dy-L*nint(dy/L)
            dz=r(i,3)-r(j,3)
            dz=dz-L*nint(dz/L)

            d2=dx**2d0+dy**2d0+dz**2d0
            d=sqrt(d2)
            if(d2.lt.cut**2)then    ! Test cutoff
                d2i=1/d2
                d6i=d2i**3

                ! ### Compute Forces ### !
                
!                ff=48d0/d**14d0-24d0/d**8d0
                 ff=(48d0*d2i*d6i*(d6i-0.5))
                force(i,1)=force(i,1)+ff*dx
                force(i,2)=force(i,2)+ff*dy
                force(i,3)=force(i,3)+ff*dz
                force(j,1)=force(j,1)-ff*dx
                force(j,2)=force(j,2)-ff*dy
                force(j,3)=force(j,3)-ff*dz
!                force(i,1)=force(i,1)+ff*dx
!                force(i,2)=force(i,2)+ff*dy
!                force(i,3)=force(i,3)+ff*dz

                ! ### Compute Potential Energy ### !
                upot=upot+4d0*(1d0/d**12d0-1d0/d**6d0)

                ! ### Compute Pressure ### !
                press=press+ff*dx*dx+ff*dy*dy+ff*dz*dz
       endif
        end if
    end do
end do
!print*,cont
call cpu_time(finishF(cont))

timeforces(cont)=finishF(cont)-startF(cont)
!media=sum(timeforces(:cont)/cont)
!print*, 'Total time=',finishF(cont)-startF(cont), 'seconds'
!print*, 'Total time=',timeforces(cont), 'seconds'
!print*, 'media', media
end subroutine forces_vlist
