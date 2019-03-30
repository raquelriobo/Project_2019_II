subroutine moment_mpi(time, vel, N, part1, part2, nini, &
root, rank)
implicit none
include 'mpif.h'
real*8 :: vel(N,3)
real*8 :: time, px, py, pz, px_aux, py_aux, pz_aux
integer :: ierr, rank, root, nini
integer :: N,i,part1,part2
open(unit=25, file="Momentum.txt")
px_aux=0
py_aux=0
pz_aux=0
px=0
py=0
pz=0

if (rank .ne. root)then
 do i=nini+1, nini+part1
  px_aux=px_aux+vel(i,1)
  py_aux=py_aux+vel(i,2)
  pz_aux=pz_aux+vel(i,3)
 end do

else

 do i=1,part2
  px_aux=px_aux+vel(i,1)
  py_aux=py_aux+vel(i,2)
  pz_aux=pz_aux+vel(i,3)
 end do
end if

call MPI_REDUCE(px_aux,px,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(py_aux,py,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)
call MPI_REDUCE(pz_aux,pz,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,ierr)

if (rank==root)then
write(25,*) time,px,py,pz
end if

end subroutine moment_mpi
