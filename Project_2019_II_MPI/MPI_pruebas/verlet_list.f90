! ### FORCE ### !
subroutine new_vlist(i,L,N,r,d2,cut,nlist,list)
implicit none
integer :: i,j,jj
integer :: N                    ! Number of particles
integer :: nlist(N-1),list(N-1,N-1)
real(8) :: L                    ! Cell long
real(8) :: r(N,3), xv(N,3)               ! Coordinates matrix
real(8) :: cut                  ! Cut-off
real(8) :: dx,dy,dz             ! Distance components
real(8) :: ff, d2, d2i, d6i
real(8) :: force(N,3)           ! Forces matrix
real(8) :: upot                 ! Potential energy
real(8) :: press                ! Pressure

print*, N

nlist=0
xv=r

do i=1,N-1
  do j=i+1,N
        dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
        dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
        dy=r(i,2)-r(j,2)
        dy=dy-L*nint(dy/L)
        dz=r(i,3)-r(j,3)
        dz=dz-L*nint(dz/L)

        d2=dx**2d0+dy**2d0+dz**2d0
        if(d2.lt.cut**2)then    ! Add to the lists
        !d2i=1/d2
        !d6i=d2i**
        nlist(i)=nlist(i)+1
        nlist(j)=nlist(j)+1
        list(i,nlist(i))=j
        list(i,nlist(j))=i
        end if
    end do
end do
end subroutine new_vlist
