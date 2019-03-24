! ### Kinetic energy ### !

subroutine kinetic_en(vel,N,kin,part1,part2,nini,nfin,nini_first,nfin_first,&
root,rank)
implicit none
include 'mpif.h'
integer :: i
integer :: N                      ! Number of particles
integer :: part1,part2
real*8 :: kin,kin_aux             ! Kinetic energy
real*8 :: vel(N,3)                ! Velocity matrix
real*8 :: vpart1(part1,3),vpart2(part2,3)
integer :: nini,nfin,nini_first,nfin_first
integer :: root,rank,ierr


kin=0.0
kin_aux=0
vpart1=0
vpart2=0

!### Separation of the velocity matrix ###!

if (rank.ne.0)then
    do i=1,part1
      vpart1(i,:)=vel(nini+i,:)
    end do
else
    do i=1,part2
      vpart2(i,:)=vel(i,:)
    end do
end if

!### Kinetic energy calculation ###!

if (rank.ne.0)then
    do i=1,N
        kin_aux=kin_aux+0.5d0*(vpart1(i,1)**2d0+vpart1(i,2)**2d0+vpart1(i,3)**2d0)
    end do
else
    do i=1,N
        kin_aux=kin_aux+0.5d0*(vpart2(i,1)**2d0+vpart2(i,2)**2d0+vpart2(i,3)**2d0)
    end do
end if

call MPI_ALLREDUCE(kin_aux,kin,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

endsubroutine kinetic_en
