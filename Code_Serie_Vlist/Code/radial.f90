!####################################!
!### Radial distribution function ###!
!####################################!

subroutine rdf(r,npart,L,switch,nhis,rho,delg,ngr,g)
implicit none
real*8 :: L                !Width of simulation box
real*8 :: r(npart,3)       !Coordinates matrix
real*8 :: rho              !Density reduced units
real*8 :: xr,yr,zr         !Distance between 2 particles components
real*8 :: delg             !Size of the bins
real*8 :: vb               !dr, volume between 2 bins
real*8 :: g(nhis)          !Radial Distribution Function
real*8 :: d                !Distance between 2 particles
real*8 :: nid              !Number of ideal gas particles in vb
real*8 :: pi
integer :: i,j
integer :: npart           !Number of particles
integer :: switch          !Switch to choose the stage of the calculation
integer :: ngr             !Counter involved in normalization of g(r)
integer :: ig              !Index of the g(r) vector
integer :: nhis            !Total number of bins, dimension of g(r)

pi=acos(-1.0d0)

!### Mode 0: initialization ###!

if (switch.eq.0)then
    ngr=0                 !Counter initialization
    delg=L/(2.0d0*nhis)   !Bin size
    g=0                   !Initialization of g(r)

!### Mode 1: Sampling during MD ###!

else if (switch.eq.1)then
    ngr=ngr+1
    do i=1,npart-1        !Loop over all pairs
       do j=i+1,npart      
           xr=r(i,1)-r(j,1)
           xr=xr-L*nint(xr/L) !PBC correction
           yr=r(i,2)-r(j,2)
           yr=yr-L*nint(yr/L)
           zr=r(i,3)-r(j,3)
           zr=zr-L*nint(zr/L)
           d=sqrt(xr**2.0d0+yr**2.0d0+zr**2.0d0)
           if (d.lt.L/2.0d0)then !Only distances within L/2 considered
               ig=int(d/delg) !Index assignment according to distance
               g(ig)=g(ig)+2  !Contribution for particle i and j
           end if
       end do
    end do

!### Mode 2: determine g(r) ###!

else if (switch.eq.2)then
   open(90,file="radial.txt")
   do i=1,nhis-1
      d=delg*(i+0.5) !Distance 
      vb=((i+1)**3.0d0-i**3.0d0)*delg**3.0d0 !Volume between bins
      nid=(4.0d0/3.0d0)*pi*vb*rho            !Number of ideal gas particles in vb
      g(i)=g(i)/(float(ngr)*npart*nid)       !Normalization of g(r)
      write(90,*) d,g(i)                     !Print
   end do
end if
end subroutine rdf 
