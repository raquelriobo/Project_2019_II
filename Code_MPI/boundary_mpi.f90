subroutine boundary_conditions(r, n_part, L, part1,part2,root,rank,nini,resizedtype,size)
implicit none
include 'mpif.h'
integer :: i
integer, intent(in) :: size
integer, intent(in) :: n_part           ! Number of particles
integer, intent(in) :: rank
integer, intent(in) :: root 
integer, intent(in) :: nini
integer, intent(in) :: resizedtype
integer, intent(in) :: part1,part2
real(8), intent(in) :: L                ! Cell longitude
real(8)             :: r(n_part,3)      ! Coordinates matrix
real*8              :: rpart1(part1,3)
real*8              :: rpart2(part2,3)
real*8              :: r_aux(part1*size,3)
integer             :: numpart
integer             :: ierr

numpart = part1*3

if(rank.ne.0) then
do i=1,part1
  rpart1(i,:) = r(nini+i,:)
end do
else
do i=1,part2
  rpart2(i,:) = r(i,:)
end do
end if

if(rank.ne.0) then
  rpart1 = rpart1-nint(rpart1/L)*L
else
  rpart2 = rpart2-nint(rpart2/L)*L
end if

call MPI_Gather(rpart1,numpart,MPI_REAL8,r_aux,1,resizedtype,root,MPI_COMM_WORLD,ierr)

if(rank.eq.0) then
do i=1,part2
  r(i,:) = rpart2(i,:)
end do
do i=1,n_part-part2
  r(i+part2,:) = r_aux(part1+i,:)
end do
end if

call MPI_Bcast(r,n_part,MPI_REAL8,root,MPI_COMM_WORLD,ierr)
end subroutine boundary_conditions
