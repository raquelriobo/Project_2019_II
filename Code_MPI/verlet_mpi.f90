! ### INTEGRATION : VERLET VELOCITY #### !
subroutine verlet_mpi(N,cut_off,press,r,v,F,dt,E_pot,L,root,rank,part1,part2,&
nini,size,resizedtype)
implicit none
include 'mpif.h'
integer, intent(in) :: N                       ! Number of part1icles
integer, intent(in) :: part1,part2                 !Number of particles bloc1, bloc2
integer, intent(in) :: size
integer, intent(in) :: resizedtype
real*8, intent(in) :: cut_off                      
real*8, intent(in) :: dt                           ! Time steps
real*8             :: r(N,3), rnew(N,3) ! Coordinates matrix
real*8             :: rpart1(part1,3),rpart2(part2,3)
real*8             :: rnew_part1(part1,3),rnew_part2(part2,3)
real*8             :: rnew_part1_aux(part1*size,3)
real*8             :: v(N,3), vnew(N,3) ! Velocity matrix
real*8             :: vpart1(part1,3),vpart2(part2,3)
real*8             :: vnew_part1(part1,3),vnew_part2(part2,3)
real*8             :: vnew_part1_aux(part1*size,3)
real*8             :: F(N,3), Fnew(N,3) ! Forces matrix
real*8             :: Fpart1(part1,3),Fpart2(part2,3)
real*8             :: Fnew_part1(part1,3),Fnew_part2(part2,3)
real*8             :: L                            ! Cell longitude
real*8             :: press                        ! Pressure
real*8             :: E_Pot                        ! Potential energy
integer             :: i,j
integer             :: ierr,rank,root
integer             :: nini
integer             :: numpart

! ### Compute new positions with verlet velocity algorithm ### !
if(rank.ne.0) then
do i=1,part1
  rpart1(i,:) = r(nini+i,:)
  vpart1(i,:) = v(nini+i,:)
  Fpart1(i,:) = F(nini+i,:)
  rnew_part1(i,:) = rpart1(i,:) + vpart1(i,:)*dt + &
                    0.5*Fpart1(i,:)*dt**2
end do
else
do i=1,part2
  rpart2(i,:) = r(i,:)
  vpart2(i,:) = v(i,:)
  Fpart2(i,:) = F(i,:)
  rnew_part2(i,:) = rpart2(i,:) + vpart2(i,:)*dt + &
                    0.5*Fpart2(i,:)*dt**2
end do
end if

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

numpart=part1*3
call MPI_Gather(rnew_part1,numpart,MPI_REAL8,rnew_part1_aux,1, &
                resizedtype,root,MPI_COMM_WORLD,ierr)

if(rank==root) then
    do i=1,part2
      rnew(i,:) = rnew_part2(i,:)
    end do
    do i=1,N-part2
      rnew(i+part2,:) = rnew_part1_aux(part1+i,:)
    end do
end if 
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ### Compute new forces ###!
CALL forces_LJ(L,N,rnew,cut_off,Fnew,press,E_pot)

! ### Compute new velocity ### !
if(rank.ne.0) then
do i=1,part1
  Fnew_part1(i,:) = Fnew(nini+i,:)
  vnew_part1(i,:) = vpart1(i,:) + 0.5*(Fnew_part1(i,:)+Fpart1(i,:))*dt
end do
else
do i=1,part2
  Fnew_part2(i,:) = Fnew(i,:)
  vnew_part2(i,:) = vpart2(i,:) + 0.5*(Fnew_part2(i,:)+Fpart2(i,:))*dt
end do
end if

call MPI_Gather(vnew_part1,numpart,MPI_REAL8,vnew_part1_aux,1,&
                resizedtype,root,MPI_COMM_WORLD,ierr)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

if(rank==root) then
    do i=1,part2
      vnew(i,:) = vnew_part2(i,:)
    end do
    do i=1,N-part2
      vnew(i+part2,:) = vnew_part1_aux(part1+i,:)
    end do
end if 
                
! ### Update #### ! 
if(rank==0) then
  r = rnew ;    v = vnew ;    F = Fnew ;
end if

call boundary_conditions(r,N,L,part1,part2,root,rank,nini,resizedtype,size) 
end subroutine verlet_mpi
