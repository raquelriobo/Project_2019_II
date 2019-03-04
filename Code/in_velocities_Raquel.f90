! ### Velocities inialization ### !
subroutine in_velocity(vel,N,Temp)
implicit none
integer :: i,j
integer :: N        ! Number of particles
real(8) :: vel(N,3) ! Velocity matrix
real(8) :: Temp     ! Temperature
real(8) :: kin      ! Kinetic Energy
real(8) :: suma(N)
real(8) :: r1

! ### Random Velocities ### !
call random_seed
do i=1,N
    call random_number(r1)
    vel(i,:)=(2d0*r1-1d0)/2d0
end do

! ### Kinetic energy ### !
call kinetic_en(vel,N,kin)

! ### Reescaling Kinetic Energy ### !
vel=vel*sqrt(N*3d0*Temp/(2d0*kin))

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
end subroutine in_velocity
