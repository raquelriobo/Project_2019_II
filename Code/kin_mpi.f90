!Subroutine to compute kinetic energy in parallel

subroutine ekin_mpi(vel,N,kintot,part1,root,rank)
implicit none
include 'mpif.h'
real*8, intent(in) :: vel(N,3)
integer :: N,i
real*8 :: kin
real*8 , intent(out) :: kintot
integer :: ierr,rank,size,root
!real*8, allocatable :: vpart1(:,:),vpart2(:,:)
integer :: part1,part2
real*8 :: vpart1(part1,N)

kintot=0
!root=0
!call MPI_INIT(ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
!call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

!if(mod(N,size).ne.0)then
!  print*, "Invalid num of processors"
!end if
  

!part1=N/(size)
!part2=N-part1

!allocate(vpart1(part1,3))
!allocate(vpart2(part2,3))

!print*, 'Process starts on processor ',rank



call MPI_SCATTER(vel,part1*3,MPI_REAL8,vpart1,part1*3,MPI_REAL8,root,MPI_COMM_WORLD,ierr)

!Calculate kin
kin=0.0d0
do i=1,part1
  kin=kin+0.5d0*(vpart1(i,1)**2.0d0+vpart1(i,2)**2.0d0+vpart1(i,3)**2.0d0)
end do

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) then
        print *, "   ==> END BARRIER TASK <== "
    end if

call MPI_ALLREDUCE(kin,kintot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

!print*,'Finish processor ',rank
!print*, kintot


!call MPI_FINALIZE(ierr)
end subroutine ekin_mpi
