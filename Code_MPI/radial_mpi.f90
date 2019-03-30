!####################################!
!### Radial distribution function ###!
!####################################!

subroutine rdf_mpi(r,npart,L,switch,nhis,rho,delg,ngr,g,rank,root,nini,part1,part2,nlist,list)
implicit none
include 'mpif.h'
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
integer :: i,j,jj
integer :: npart           !Number of particles
integer :: switch          !Switch to choose the stage of the calculation
integer :: ngr             !Counter involved in normalization of g(r)
integer :: ig              !Index of the g(r) vector
integer :: nhis            !Total number of bins, dimension of g(r)

integer :: ierr,rank,root
integer :: nini
integer :: part1, part2
integer :: g_part(nhis),g_aux(nhis)
integer :: nlist(npart),list(npart,npart-1)

pi=acos(-1.0d0)

!### Mode 0: initialization ###!

if (switch.eq.0)then
    ngr=0                 !Counter initialization
    delg=(L*0.4)/nhis   !Bin size
    g=0                   !Initialization of g(r)

!### Mode 1: Sampling during MD ###!

else if (switch.eq.1)then
g_part=0
if (rank.ne.root) then
    do i=nini+i,nini+part1        
       do jj=1,nlist(i)
           j=list(i,jj)      
           xr=r(i,1)-r(j,1)
           xr=xr-L*nint(xr/L) !PBC correction
           yr=r(i,2)-r(j,2)
           yr=yr-L*nint(yr/L)
           zr=r(i,3)-r(j,3)
           zr=zr-L*nint(zr/L)
           d=sqrt(xr**2.0d0+yr**2.0d0+zr**2.0d0)
           if (d.lt.L*0.4)then !Only distances within L/2 considered
               ig=int(d/delg) !Index assignment according to distance
               g_part(ig)=g_part(ig)+1  !Contribution for particle i
           end if
       end do
    end do
else
    ngr=ngr+1
    do i=1,part2     
       do jj=1,nlist(i)
           j=list(i,jj)
           xr=r(i,1)-r(j,1)
           xr=xr-L*nint(xr/L) !PBC correction
           yr=r(i,2)-r(j,2)
           yr=yr-L*nint(yr/L)
           zr=r(i,3)-r(j,3)
           zr=zr-L*nint(zr/L)
           d=sqrt(xr**2.0d0+yr**2.0d0+zr**2.0d0)
           if (d.lt.L*0.4)then !Only distances within L/2 considered
               ig=int(d/delg) !Index assignment according to distance
               g_part(ig)=g_part(ig)+1  !Contribution for particle i
           end if
       end do
    end do
end if

call MPI_REDUCE(g_part,g_aux,nhis,MPI_INTEGER,MPI_SUM,root,MPI_COMM_WORLD,&
ierr)
if(rank.eq.root)then
   g=g+g_aux
end if

!### Mode 2: determine g(r) ###!

else if (switch.eq.2.and.rank.eq.root)then
   open(90,file="radial.txt")
   do i=1,nhis-1
      d=delg*(i+0.5) !Distance 
      vb=((i+1)**3.0d0-i**3.0d0)*delg**3.0d0 !Volume between bins
      nid=(4.0d0/3.0d0)*pi*vb*rho            !Number of ideal gas particles in vb
      g(i)=g(i)/(float(ngr)*npart*nid)       !Normalization of g(r)
      write(90,*) d,g(i)                     !Print
   end do
end if
end subroutine rdf_mpi 
