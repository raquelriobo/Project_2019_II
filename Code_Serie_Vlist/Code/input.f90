subroutine input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime)
implicit none
real*8, parameter :: kb = 1.3806d-23    !Boltzmann constant            
integer           :: M                  !Number particles/unit cell
integer           :: nhis               !Number steps to compute g(r)
real*8            :: dt                 !Time step
real*8            :: mass               
real*8            :: density
real*8            :: Temp               !Temperature
real*8            :: Maxtime            !Maximum time reached 
real*8            :: press              !Pressure
real*8            :: sigma, eps         !L-J parameters
real*8            :: Na = 6.022d23      !NAvogadro        
read(5,*)
read(5,*)
read(5,*) M, dt, mass, Temp, press, sigma, eps, nhis, Maxtime
 
press = press/9.869d-6                      !Convertion from atm to Pa
density = press/(kb*Temp)                   !Computation of density  

! ### Transformation to reduced units ### !
density = density * sigma**3d0              !Transformation to reduced units
Temp = (kb*Temp*Na)/eps                     !Transformation to reduced units

end subroutine input 
