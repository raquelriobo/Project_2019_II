
subroutine units_print(time,potential_E,kinetic_E,pressure,L,dt,sigma,epsil,density,Temp)
implicit none
real*8, intent(in) :: time,potential_E,kinetic_E,pressure,dt,L,Temp
real*8 :: sigma,epsil,density
integer :: i
real*8 :: ekin, epot, pressu, p !les noves variables amb el canvi d'unitats
real*8, parameter ::Avognum=6.022140875d23

open(unit=24, file="Results.txt")


!Falta cambio de unidades de tiempo y de temperatura!!!!!
ekin=kinetic_E*epsil*1d-3
epot=potential_E*epsil*1d-3
pressu=pressure/2d0
p=density*Temp*epsil*(sigma**(-3d0))*Avognum**(-1d0)+((1d0/(3d0*(L*sigma)**3d0))&
&*epsil*Avognum**-(1d0))*pressu
 
write(24,*) time, epot, ekin, (epot+ekin), p

end subroutine units_print
