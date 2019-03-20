! ### INTEGRATION : VERLET VELOCITY #### !
subroutine verlet_mpi(n_part,cut_off,press,r,v,F,dt,E_pot,L,root,part1,rank)
implicit none
include 'mpif.h'
integer, intent(in) :: n_part                       ! Number of part1icles
real(8), intent(in) :: cut_off                      
real(8), intent(in) :: dt                           ! Time steps
real(8)             :: r(n_part,3), r_new(n_part,3) ! Coordinates matrix 
real(8)             :: v(n_part,3), v_new(n_part,3) ! Velocity matrix
real(8)             :: F(n_part,3), F_new(n_part,3) ! Forces matrix
real(8)             :: L                            ! Cell longitude
real(8)             :: press                        ! Pressure
real(8)             :: E_Pot                        ! Potential energy
integer             :: i,j
integer             :: part1
real(8)             :: vpart1(part1,3),fpart1(part1,3),rpart1(part1,3)
integer             :: ierr,rank,root
real(8)             :: r_new_part(n_part,3)
integer             :: dimens
integer             :: numworkers !Número de procesadores
integer             :: oldsize(2)
integer             :: newsize(2)
integer             :: block
integer             :: starts(2)
integer             :: dtype(n_part*3)
integer             :: nblock_col,nblock_fil

rpart1=0
vpart1=0
fpart1=0
r_new=0

numworkers=4 !Lo tendría que pasar el programa desde el main
dimens=n_part*3 !Dimensión de las matrices originales
oldsize(1) = n_part  !Filas originales
oldsize(2) = 3  !Columnas originales
newsize(1) = int(n_part/numworkers) !Filas parciales
newsize(2) = 3 !Columnas parciales
nblock_fil=numworkers
nblock_col=1

block=0

do i=1,nblock_fil
        do j=1,nblock_col
           ! Creating data type subblocs
            starts(1) = (i-1)*newsize(1)
            starts(2) = (j-1)*newsize(2)
            call MPI_TYPE_CREATE_SUBARRAY(2,oldsize,newsize,starts,MPI_ORDER_FORTRAN, &
                                          MPI_REAL8,dtype(block),ierr) !create subarray
            call MPI_TYPE_COMMIT(dtype(block),ierr) !commmits a data type

            block = block + 1
        end do
end do

print*, 'Process starts on processor ',rank


! Scatter old velocities
!call MPI_SCATTERv(v,n_part*3,part1,MPI_REAL8,vpart1,part1*3,MPI_REAL8,root,MPI_COMM_WORLD,ierr)

! Scatter old positions
!call MPI_SCATTERv(r,n_part*3,part1,MPI_REAL8,rpart1,part1*3,MPI_REAL8,root,MPI_COMM_WORLD,ierr)

! Scatter old forces
!call MPI_SCATTERv(f,n_part*3,part1,MPI_REAL8,fpart1,part1*3,MPI_REAL8,root,MPI_COMM_WORLD,ierr)

! Barrier so all processors have all info
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) then
        print *, "   ==> END BARRIER TASK <== "
    end if
              
! ### Compute new positions with verlet velocity algorithm ###!
!Each processor: 
do i=1,part1
  do j=1,3
    r_new_part(i,j) = rpart1(i,j) + vpart1(i,j)*dt + &
    0.5*Fpart1(i,j)*dt**2
  end do
end do 

print*,r_new_part(1,1),r_new_part(1,2),r_new_part(1,3)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) then
        print *, "   ==> END BARRIER TASK <== "
    end if

call MPI_ALLGATHERv(r_new_part,part1*3,MPI_REAL8,r_new,n_part*3,part1,MPI_REAL8,MPI_COMM_WORLD,ierr)
!call MPI_ALLREDUCE(r_new_part,r_new,n_part*3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

do i=1,n_part
  print*, r_new(i,1),r_new(i,2),r_new(i,3)
end do


call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank .eq. 0) then
        print *, "   ==> END BARRIER TASK <== "
    end if


print*, 'Process finishes on processor ',rank

!print*, r_new(1,1),r_new(1,2),r_new(1,3)

! ### Compute new forces ###!
CALL forces_LJ(L,n_part,r_new,cut_off,F_new,press,E_pot)

! ### Compute new velocity ###! 
do i=1,n_part
    v_new(i,:) = v(i,:) + &
    0.5*(F(i,:)+F_new(i,:))*dt
end do
                
! ### Update #### ! 
r = r_new ;    v = v_new ;    F = F_new ;
r_new = 0.;    v_new = 0.;    F_new = 0.;

call boundary_conditions(r,n_part,L) 
end subroutine verlet_mpi
