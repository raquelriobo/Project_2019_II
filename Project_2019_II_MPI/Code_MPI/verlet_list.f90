! ### FORCE ### !
subroutine new_vlist(i,L,N,r,cut,nlist,list,rold,root,rank,part1,part2,nini,size)
implicit none
include 'mpif.h' 
integer :: i,j,jj
integer :: N                    ! Number of particles
integer :: nlist(N),list(N,N-1) !!!!!!!!******** Seria list(N,N-1) PERO NO VA
real(8) :: L                    ! Cell long
real(8) :: r(N,3), rold(N,3)               ! Coordinates matrix
real(8) :: cut                  ! Cut-off
real(8) :: dx,dy,dz             ! Distance components
real(8) :: ff, d2, d2i, d6i
real(8) :: force(N,3)           ! Forces matrix
real(8) :: upot                 ! Potential energy
real(8) :: press                ! Pressure
integer :: ierr,rank,root
integer :: nini,size
integer :: part1, part2
integer :: list_aux(part1*size,N-1), nlist_aux(part1*size)
integer :: nlist_1(part1),nlist_2(part2)
integer :: list_1(part1,N-1),list_2(part2,N-1)
integer :: sizes(2),subsizes(2),starts(2)
integer :: blocktype,resizedtype
integer :: intsize
integer (kind=MPI_Address_kind) :: start, extent

list_1=0
list_2=0
nlist_1=0
nlist_2=0
rold=r
nlist_aux=0
list_aux=0

if (rank.ne.root) then
    do i=nini+1,nini+part1
        do j=1,N
            if(i.ne.j)then
                dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
                dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
                dy=r(i,2)-r(j,2)
                dy=dy-L*nint(dy/L)
                dz=r(i,3)-r(j,3)
                dz=dz-L*nint(dz/L)

                d2=dx**2d0+dy**2d0+dz**2d0
                    if(d2.lt.cut**2)then    ! Add to the lists
                        nlist_1(i-nini)=nlist_1(i-nini)+1
                        list_1(i-nini,nlist_1(i-nini))=j
                    end if
            endif
        end do
    end do
else
    do i=1,part2
        do j=1,N
            if(i.ne.j)then
                dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
                dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
                dy=r(i,2)-r(j,2)
                dy=dy-L*nint(dy/L)
                dz=r(i,3)-r(j,3)
                dz=dz-L*nint(dz/L)

                d2=dx**2d0+dy**2d0+dz**2d0
                    if(d2.lt.cut**2)then    ! Add to the lists
                        nlist_2(i)=nlist_2(i)+1
                        list_2(i,nlist_2(i))=j
                    end if
            endif
        end do
    end do
end if

!if(rank.eq.1)then
!do i=1,nlist_1(108)
!print*,nlist_1(1)
!  print*,list_1(108,:)
!end do
!end if

sizes=  [ part1*size, N-1 ] !sizes array total
subsizes = [ part1, N-1 ] !subsizes subarrays
starts   = [ 0, 0 ]


call MPI_Type_create_subarray( 2, sizes, subsizes, starts,     &
                                 MPI_ORDER_FORTRAN, MPI_INTEGER, &
                                 blocktype, ierr)
start = 0
call MPI_Type_size(MPI_INTEGER, intsize, ierr)
extent = intsize * part1

call MPI_Type_create_resized(blocktype, start, extent, resizedtype, ierr)
call MPI_Type_commit(resizedtype,ierr)

do i=1,nlist_1(108)
   print*, list_1(108,i)
end do
call MPI_BARRIER(MPI_COMM_WORLD,ierr)


call MPI_Gather(list_1,part1*(N-1),MPI_INTEGER,list_aux,1,MPI_INTEGER,root,MPI_COMM_WORLD, ierr)

!call MPI_Gather(nlist_1,1,MPI_INTEGER,nlist_aux,1,MPI_INTEGER,root,&
!                  MPI_COMM_WORLD, ierr)

if (rank==root)then
 do i =1,part2
   list(i,:)=list_2(i,:)
 end do
 do i=1,N-part2
    list(i+part2,:)=list_aux(part1+i,:)
 end do
 do i =1,part2
   nlist(i)=nlist_2(i)
 end do
 do i=1,N-part2
    nlist(i+part2)=nlist_aux(part1+i)
 end do
end if

call MPI_BCAST(list,N*(N-1), MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nlist, N, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)

end subroutine new_vlist
