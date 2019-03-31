! ### Total Momentum MPI ### !
subroutine moment_mpi(time, vel, N, part1, part2, nini, &
root, rank)
implicit none
include 'mpif.h'
real(8) :: vel(N,3)                ! Velocity Matrix
real(8) :: time                    ! Time 
real(8) :: px, py, pz              ! Momentum components 
real(8) :: px_aux, py_aux, pz_aux  ! Auxiliar momentum components
integer :: ierr, rank, root 
integer :: nini                    ! Index MPI data partition
integer :: N                       ! Number of particles
integer :: i
integer :: part1                   ! Number of particles to send to processor
integer :: part2                   ! Number of particles to send to root

px_aux=0
py_aux=0
pz_aux=0
px=0
py=0
pz=0

if (rank .ne. root)then

! Send data to processors !
 do i=nini+1, nini+part1
  px_aux=px_aux+vel(i,1)
  py_aux=py_aux+vel(i,2)
  pz_aux=pz_aux+vel(i,3)
 end do

else

! Send data to root !
 do i=1,part2
  px_aux=px_aux+vel(i,1)
  py_aux=py_aux+vel(i,2)
  pz_aux=pz_aux+vel(i,3)
 end do
end if

! Add all terms for every momentum component !
call MPI_REDUCE(px_aux,px,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(py_aux,py,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(pz_aux,pz,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)

! ### Print results ### !
if (rank==root)then
open(unit=27, file="Momentum_xy.txt")
open(unit=28, file="Momentum_z.txt")

write(27,*) time,px,py
write(28,*) time,pz
end if

end subroutine moment_mpi
