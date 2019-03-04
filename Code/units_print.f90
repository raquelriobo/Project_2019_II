! ### Subroutine that change reduced units and print the results ### !
subroutine units_print(time,potential_E,kinetic_E,pressure,L,dt,sigma,epsil,density,Temp,mass)
implicit none
integer             :: i
real(8), intent(in) :: density, L, sigma,epsil
real(8), intent(in) :: time, dt
real(8), intent(in) :: potential_E,kinetic_E ! Energy
real(8), intent(in) :: pressure,Temp,mass
real(8)             :: time2, Temp2
real(8)             :: ekin, epot, pressure2 ! New variables with non-reduced units
real(8), parameter  :: Avognum = 6.022140875d23
real(8), parameter  :: kb = 1.38064852d-23

open(unit=24, file="Results.txt")

ekin=kinetic_E*epsil*1d-3               ! Kinetic Energy
epot=potential_E*epsil*1d-3             ! Potential Energy

pressure2=pressure/2d0
pressure2=density*Temp*epsil*(sigma**(-3d0))*Avognum**(-1d0)+((1d0/(3d0*(L*sigma)**3d0))&
&*epsil*Avognum**(-1d0))*pressure2     ! Pressure

time2=time/(sigma*(sqrt((mass*(1d-3))/epsil)))
time2=time2*(1d-12)                     ! Time in picoseconds

Temp2=(Temp*epsil)/(kb*Avognum)         ! Temperature

! Print the results !
write(24,*) time2, epot, ekin, (epot+ekin), pressure2, Temp2
! time, Potential Energy, Kinetic Energy, Total Energy, Pressure and Temperature !

end subroutine units_print
