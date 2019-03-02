! ********** FORCE ********** !
subroutine forces_LJ(L,N,r,cut,upot,force)
implicit none
integer :: i,j,N
real*8 :: L,r(N,3),upot,dx,dy,dz,d2,cut,force(N,3),d, d2i,d6i,ff

upot=0d0
force=0d0

! **** Periodic boundary conditions **** !
do i=1,N-1
  do j=1+i,N 			! * Loop over all pairs
!    if(i.ne.j)then
      dx=r(i,1)-r(j,1)   	! * MINUMUM IMAGE CONVENTION donde estaria la vecina mas proxima
      dx=dx-L*nint(dx/L)	! * PBC
      dy=r(i,2)-r(j,2)
      dy=dy-L*nint(dy/L)
      dz=r(i,3)-r(j,3)
      dz=dz-L*nint(dz/L)
      d2=dx**2d0+dy**2d0+dz**2d0
     ! d=sqrt(d2)		! * Distance between particles using PBC
      if(d2.lt.cut**2)then	! * Test cutoff
	d2i=1/d2
	d6i=d2i**3
!	ff=(48d0/d**14d0-24d0/d**8d0)
	ff=(48d0*d2i*d6i*(d6i-0.5))
        force(i,1)=force(i,1)+ff*dx
        force(i,2)=force(i,2)+ff*dy
        force(i,3)=force(i,3)+ff*dz
        force(j,1)=force(j,1)-ff*dx
        force(j,2)=force(j,2)-ff*dy
        force(j,3)=force(j,3)-ff*dz
	
        !upot=upot+4d0*(1d0/d**12d0-1d0/d**6d0)
        upot=upot+4d0*d6i*(d6i-1)
!      end if
    end if
  end do
end do
end subroutine forces_LJ
