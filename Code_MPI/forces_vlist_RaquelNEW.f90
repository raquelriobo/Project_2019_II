! ### FORCE ### !
subroutine forces_vlist(L,N,r,cut,force,press,upot,nlist,list,rank,root,part1,part2,nini,size,resizedtype)

implicit none
include 'mpif.h'
integer :: i,j,jj
integer :: N                    ! Number of particles
integer :: nlist(N),list(N,N-1)
real(8) :: L                    ! Cell long
real(8) :: r(N,3)               ! Coordinates matrix
real(8) :: cut                  ! Cut-off
real(8) :: dx,dy,dz             ! Distance components
real(8) :: ff, d2, d2i, d6i, d
real(8) :: force(N,3)           ! Forces matrix
real(8) :: upot                 ! Potential energy
real(8) :: press                ! Pressure
integer :: size,ierr,rank,root
integer :: nini,resizedtype
integer :: part1, part2
real*8  :: fpart1(part1,3),fpart2(part2,3)
real*8 :: press_aux,upot_aux,f_aux(part1*size,3)

force=0d0
upot=0d0
upot_aux=0d0
press=0d0
press_aux=0d0
fpart1=0
fpart2=0
f_aux=0

!### Partial matrixes ###!
!if (rank.ne.root) then
!    do i=1,part1
!        fpart1(i,:)=force(nini+i,:)
!        print*,fpart1(i,:)
!    end do
!else
!    do i=1,part2
!        fpart2(i,:)=force(i,:)
!    end do
!end if
    


if (rank.ne.root) then
    do i=nini+1,nini+part1
        do jj=1,nlist(i)
            j=list(i,jj)
            dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
            dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
            dy=r(i,2)-r(j,2)
            dy=dy-L*nint(dy/L)
            dz=r(i,3)-r(j,3)
            dz=dz-L*nint(dz/L)

            d2=dx**2d0+dy**2d0+dz**2d0
            d=sqrt(d2)
            if(d2.lt.cut**2)then    ! Test cutoff
                d2i=1/d2
                d6i=d2i**3

                ! ### Compute Forces ### !
                
                ff=48d0/d**14d0-24d0/d**8d0
                
                fpart1(i-nini,1)=fpart1(i-nini,1)+ff*dx
                fpart1(i-nini,2)=fpart1(i-nini,2)+ff*dy
                fpart1(i-nini,3)=fpart1(i-nini,3)+ff*dz

                ! ### Compute Potential Energy ### !
                upot_aux=upot_aux+2d0*(1d0/d**12d0-1d0/d**6d0)

                ! ### Compute Pressure ### !
                press_aux=press_aux+ff*dx*dx+ff*dy*dy+ff*dz*dz
        end if
    end do
end do
else
    do i=1,part2
        do jj=1,nlist(i)
            j=list(i,jj)
            dx=r(i,1)-r(j,1)        ! Minimum Image Convention (where the cosest neighbor is)
            dx=dx-L*nint(dx/L)      ! Periodic Boundary Conditions
            dy=r(i,2)-r(j,2)
            dy=dy-L*nint(dy/L)
            dz=r(i,3)-r(j,3)
            dz=dz-L*nint(dz/L)

            d2=dx**2d0+dy**2d0+dz**2d0
            d=sqrt(d2)
            if(d2.lt.cut**2)then    ! Test cutoff
                d2i=1/d2
                d6i=d2i**3

                ! ### Compute Forces ### !

                ff=48d0/d**14d0-24d0/d**8d0

                fpart2(i,1)=fpart2(i,1)+ff*dx
                fpart2(i,2)=fpart2(i,2)+ff*dy
                fpart2(i,3)=fpart2(i,3)+ff*dz

                ! ### Compute Potential Energy ### !
                upot_aux=upot_aux+2d0*(1d0/d**12d0-1d0/d**6d0)

                ! ### Compute Pressure ### !
                press_aux=press_aux+ff*dx*dx+ff*dy*dy+ff*dz*dz
        end if
    end do
end do
end if

call MPI_Gather(fpart1,part1*3,MPI_REAL8,f_aux,1,resizedtype,root,&
MPI_COMM_WORLD,ierr)

!if (rank.eq.root)then
!print*, f_aux(67,1),f_aux(67,2),f_aux(67,3)
!end if

call MPI_ALLREDUCE(press_aux,press,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_ALLREDUCE(upot_aux,upot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

if (rank.eq.root)then
 do i =1,part2
   force(i,:)=fpart2(i,:)
 end do
 do i=1,N-part2
    force(i+part2,:)=f_aux(part1+i,:)
 end do
end if

!if (rank.eq.root)then
!print*, force(107,1),force(107,2),force(107,3)
!end if

call MPI_BCAST(force, N*3, MPI_REAL8, root, MPI_COMM_WORLD, ierr)

!print*, force(107,1),force(107,2),force(107,3)


end subroutine forces_vlist
