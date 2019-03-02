!Radial distribution function
      subroutine rdf(r,npart,L,switch,nhis,rho,delg,ngr,g)
      implicit none
      real*8 :: L,r(npart,3),rho
      real*8 :: xr,yr,zr,delg,vb,g(nhis),d,pi,nid
      integer :: i,j,npart,switch,ngr,ig,nhis
      pi=acos(-1.0d0)
      if (switch.eq.0)then
          ngr=0
          delg=L/(2.0d0*nhis)
          g=0
      else if (switch.eq.1)then
          ngr=ngr+1
          do i=1,npart-1
             do j=i+1,npart
                 xr=r(i,1)-r(j,1)
                 yr=r(i,2)-r(j,2)
                 zr=r(i,3)-r(j,3)
                 d=sqrt(xr**2.0d0+yr**2.0d0+zr**2.0d0)
                 if (d.lt.L/2.0d0)then
                     ig=int(d/delg)
                     g(ig)=g(ig)+2
                 end if
             end do
          end do
      else if (switch.eq.2)then
         open(90,file="radial.txt")
         do i=1,nhis
            d=delg*(i+0.5)
            vb=((i+1)**3.0d0-i**3.0d0)*delg**3.0d0
            nid=(4.0d0/3.0d0)*pi*vb*rho
            g(i)=g(i)/(ngr*npart*nid)
            write(90,*) d,g(i)
         end do
      end if
      end subroutine rdf 
