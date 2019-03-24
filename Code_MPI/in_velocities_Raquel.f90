! ### Velocities inialization ### !
subroutine in_velocity(vel,N,Temp,part1,part2,nini,nfin,nini_first,nfin_first,&
root,rank)
implicit none
integer :: i,j
integer :: N        ! Number of particles
real(8) :: vel(N,3) ! Velocity matrix
real(8) :: Temp     ! Temperature
real(8) :: kin      ! Kinetic Energy
real(8) :: suma(3)
real(8) :: r1
integer :: part1,part2,nini,nfin,nini_first,nfin_first,root,rank

! ### Random Velocities ### !

!Assign a random velocity component between -0.5 and 0.5
call random_seed
do i=1,N
    do j=1,3
      call random_number(r1)
      vel(i,j)=(2d0*r1-1d0)/2d0
    end do
end do

! ### Total velocity = 0 ### !
suma=0d0
do i=1,3
    do j=1,N
        suma(i)=suma(i)+vel(j,i)
    end do
enddo

suma=suma/N

do i=1,N
    do j=1,3
        vel(i,j)=vel(i,j)-suma(j)
    enddo
enddo

!Kinetic energy calculation
call kinetic_en(vel,N,kin,part1,part2,nini,nfin,nini_first,nfin_first,&
root,rank)

!Reescaling velocities to match kinetic energy
vel=vel*sqrt(N*3d0*Temp/(2d0*kin))

end subroutine in_velocity
