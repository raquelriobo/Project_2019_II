subroutine checkvlist(i,L,N,r,cut,cut2,nlist,list,rold,dx,dy,dz,counter,cont)
implicit none
integer :: N
integer :: i,counter,cont
integer :: nlist(N),list(N,N-1)
real*8  :: cut            !LJ forces calculation distance cut-off
real*8  :: cut2
real*8  :: L
real*8  :: resta
real*8  :: r(N,3)
real*8 :: rold(N,3)
real(8) :: dx,dy,dz

        do i=1,N
            dx=r(i,1)-rold(i,1)
            dx=dx-L*nint(dx/L)
            dy=r(i,2)-rold(i,2)
            dy=dy-L*nint(dy/L)
            dz=r(i,3)-rold(i,3)
            dz=dz-L*nint(dz/L)
            resta=dx**2d0+dy**2d0+dz**2d0
            if (resta.ge.((cut2-cut)**2d0)) then
                cont=cont+1
!                print*,counter,cont
                call new_vlist(i,L,N,r,cut2,nlist,list,rold)
                exit
            endif
        enddo
endsubroutine checkvlist
