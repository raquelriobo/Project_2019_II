! ### Velocities inicialization MPI ### !
subroutine in_velocity_mpi(vel,N,Temp,part1,part2,nini,nfin,&
root,rank,resizedtype,size)
implicit none
include 'mpif.h'
integer :: size
integer :: part1                    ! Number of particles send to procesor
integer :: part2                    ! Number of particles send to root
real(8) :: vel(N,3)                 ! Velocity Matrix
real(8) :: vel_final(part1*size,3)  ! Final Velocity Matrix 
real(8) :: veltot(part1*size,3)     ! Velocity Matrix output first gatherMPI
real(8) :: vel_end(N,3)             ! Velocity Matrix with total velocity = 0
integer :: N                        ! Number of particles
integer :: i,j
integer :: nini,nfin                ! Index MPI data partition
integer :: nini_first, nfin_first   ! Index MPI data partition
real(8) :: kin                      ! Kinetic Energy
real(8) :: r1                       ! Random number
integer :: ierr,rank,root
real(8) :: vpart1(part1,3)          ! part of velocity matrix (processor) 
real(8) :: vpart2(part2,3)          ! part of velocity matrix (root)
integer :: resizedtype
real(8) :: suma(3) 
real(8) :: Temp                     ! Temperature

! ### Random Velocities ### !
!Assign a random velocity component between -0.5 and 0.5 !

call random_seed
vpart1=0
vpart2=0

if (rank .ne. root)then

 ! Create submatrix "vpart1" !
 ! Data to send to processors !
 do i=1,part1
  do j=1,3
   call random_number(r1)
   vpart1 (i,j)=(2d0*r1-1d0)/2d0
  end do
 end do

else 

 ! Create submatrix "vpart2" !
 ! Data to send to root !
 do i=1,part2
  do j=1,3
   call random_number(r1)
   vpart2 (i,j)=(2d0*r1-1d0)/2d0
  end do
 end do

end if


call MPI_Gather( vpart1, part1*3, MPI_REAL8, &  
                     veltot, 1, resizedtype,root, &
                   MPI_COMM_WORLD, ierr)

! Join the velocity matrix !
if (rank==root)then
 do i =1,part2 
   vel_end(i,:)=vpart2(i,:)
 end do
 do i=1,N-part2
    vel_end(i+part2,:)=veltot(part1+i,:)
 end do

! ### Total velocity = 0 ### !!!!
suma=0d0
do i=1,3
    do j=1,N
        suma(i)=suma(i)+vel_end(j,i)
    end do
end do

suma=suma/N

do i=1,N
    do j=1,3
       vel_end(i,j)=vel_end(i,j)-suma(j)
    end do
end do

end if

call MPI_BCAST(vel_end, N*3, MPI_REAL8, root, MPI_COMM_WORLD, ierr)

! ### Kinetic Energy Calculation ### !
! use vel_end as input, and kin as output !
call kinetic_en(vel_end,N,kin,part1,part2,nini,nfin,nini_first,nfin_first,&
root,rank)

! ### Reescaling velocities to match kinetic energy ### !
vpart1=0
vpart2=0

if (rank .ne. root) then

 ! Create submatrix "vpart1" !
 ! Data to send to processors !
 do i=1,part1
  do j=1,3
   vpart1(i,j)=vel_end(nini+i,j)
   vpart1(i,j)=vpart1(i,j)*sqrt(N*3d0*Temp/(2d0*kin))
  end do
 end do

else

 ! Create submatrix "vpart2" !
 ! Data to send to root !
 do i=1,part2
  do j=1,3
   vpart2(i,j)=vel_end(i,j)
   vpart2(i,j)=vpart2(i,j)*sqrt(N*3d0*Temp/(2d0*kin))
   end do
 end do

end if

 call MPI_Gather(vpart1, part1*3, MPI_REAL8, &
                 vel_final, 1, resizedtype, root,&
                 MPI_COMM_WORLD, ierr)

! Join the velocity matrix !
if (rank==root)then;

 do i =1,part2 
   vel(i,:)=vpart2(i,:)
 end do

 do i=1,N-part2
    vel(i+part2,:)=vel_final(part1+i,:)
 end do

end if

! Send Velocity Matrix data !
call MPI_BCAST(vel, N*3, MPI_REAL8, root, MPI_COMM_WORLD, ierr) 


end subroutine in_velocity_mpi
