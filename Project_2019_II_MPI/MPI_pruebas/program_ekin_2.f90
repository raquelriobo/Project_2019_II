program main
implicit none
include 'mpif.h'
integer,parameter :: N=4
real*8 :: vel(N,3)
real*8 :: kin
integer :: i,ierr,rank,size,root,part1

root=0
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

part1=N/(size) !numero de particulas /num procesadores

vel(1,:)=1
vel(2,:)=2
vel(3,:)=3
vel(4,:)=4

call ekin_mpi(vel,N,kin,part1,root,rank,size)

write(6,*) "kin=",kin

end program main

!Subroutine to compute kinetic energy in parallel


subroutine ekin_mpi(vel,N,kintot,part1,root,rank,size)
implicit none
include 'mpif.h'
real*8, intent(in) :: vel(N,3)
integer :: N,i
real*8 :: kin
real*8 , intent(out) :: kintot
integer :: ierr,rank,size,root
integer :: part1
real*8  :: vpart1(part1,3)
integer :: sizes(2),subsizes(2),starts(2)
integer :: blocktype,resizedtype
integer :: intsize
integer (kind=MPI_Address_kind) :: start, extent

kintot=0

call MPI_BCAST(vel,N*3,MPI_REAL8,root,MPI_COMM_WORLD,ierr)

write(6,*) "Task",rank,vel(1,1),vel(1,2),vel(1,3)
write(6,*) "Task",rank,vel(2,1),vel(2,2),vel(2,3)
write(6,*) "Task",rank,vel(3,1),vel(3,2),vel(3,3)
write(6,*) "Task",rank,vel(4,1),vel(4,2),vel(4,3)

!Calculate kin
kin=0.0d0
do i=rank*part1+1,rank*part1+part1 ! numero del procesador(task ID) 0*4part/2procesadores+2=1 a 0*1+2=2
  kin=kin+0.5d0*(vel(i,1)**2.0d0+vel(i,2)**2.0d0+vel(i,3)**2.0d0)
  print*,rank,kin
end do

call MPI_ALLREDUCE(kin,kintot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

call MPI_FINALIZE(ierr)
end subroutine ekin_mpi
