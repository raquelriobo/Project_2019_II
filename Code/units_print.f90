! ### Subroutine that change reduced units and print the results ### !
subroutine units_print(time,potential_E,kinetic_E,pressure,L,dt,sigma,epsil,density,Temp,mass)
implicit none
real*8, intent(in) :: time,potential_E,kinetic_E,pressure,dt,L,Temp,mass
real*8 :: sigma,epsil,density, timee, Tempp, timefs
integer :: i
real*8 :: ekin, epot, pressu, p ! New variables with non-reduced units
real*8, parameter ::Avognum=6.022140875d23
real*8, parameter :: kb=1.38064852d-23
open(unit=24, file="Results.txt")

ekin=kinetic_E*epsil*1d-3 !Kinetic Energy
epot=potential_E*epsil*1d-3 !Potential Energy

pressu=pressure/2d0
p=density*Temp*epsil*(sigma**(-3d0))*Avognum**(-1d0)+((1d0/(3d0*(L*sigma)**3d0))&
&*epsil*Avognum**(-1d0))*pressu  !Pressure

timee=time/(sigma*(sqrt((mass*(1d-3))/epsil)))
timefs=timee*(1d-12)  !Time in picoseconds

Tempp=(Temp*epsil)/(kb*Avognum) !Temperature

!Print the results
write(24,*) timefs, epot, ekin, (epot+ekin), p, Tempp
! time, Potential Energy, Kinetic Energy, Total Energy, Pressure and Temperature !

end subroutine units_print
