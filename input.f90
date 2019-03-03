    subroutine input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime)
        implicit none
        real*8, parameter :: kb = 1.3806d-23
        integer           :: M, nhis
        real*8            :: dt, mass, density, Temp, Maxtime
        real*8            :: press       !Initial pression
        real*8            :: sigma
        real*8            :: eps
        
        read(5,*)
        read(5,*)
        read(5,*) M, dt, mass, Temp, press, sigma, eps, nhis, Maxtime
 

        density = (press/9.869d-6)/(kb*Temp)
        density = density * sigma**3d0 !Reduced units
        Temp = (kb*Temp*6.022d23)/eps  !Reduced units
    end subroutine input 
