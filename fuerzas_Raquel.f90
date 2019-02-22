! ********** FORCE ********** !
subroutine forces_LJ_Press(L,N,r,cut,upot,force,press)
implicit none
integer :: i,j,N
real*8 :: L,r(N,3),upot,dx,dy,dz,d2,d6,d12,cut,force(N,3),d,press,fnew

upot=0d0
force=0d0
press=0d0


do i=1,N
  do j=1,N !Loop over all pairs
    if(i.ne.j)then
      dx=r(i,1)-r(j,1)
      dx=dx-nint(dx/L)*L !Correction of distance according to PBC
      dy=r(i,2)-r(j,2)
      dy=dy-nint(dy/L)*L
      dz=r(i,3)-r(j,3)
      dz=dz-nint(dz/L)*L
      d2=dx*dx+dy*dy+dz*dz
      d=sqrt(d2) !Distance between particles using PBC
      if(sqrt(d2).lt.cut)then
	fnew=48d0/d**14d0-24d0/d**8d0 !Force actualization
	force(i,1)=force(i,1)+fnew*dx
	force(i,2)=force(i,2)+fnew*dy
	force(i,3)=force(i,3)+fnew*dz

	press=press+fnew*dx*dx+fnew*dy*dy+fnew*dz*dz

	upot=upot+2d0*(1d0/d**12d0-1d0/d**6d0)
      end if
    end if
  end do
end do
end subroutine forces_LJ_Press
