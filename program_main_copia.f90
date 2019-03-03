!Programa Principal


program main
implicit none
integer :: M,N,i,nhis,ngr,counter
real*8 :: upot,time,density,cut,L,dt,Temp,press,kin,delg
real*8 :: mass,eps,sigma,Maxtime
real*8,allocatable :: r(:,:),vel(:,:),force(:,:),g(:)


call input(M,dt,mass,density,Temp,sigma,eps,nhis,Maxtime)

N=M**(3)*4 !Number of particles
allocate(r(N,3),vel(N,3),force(N,3),g(nhis))


!System initialization

call coordenadas(N,M,r,density,L)
call boundary_conditions(r,N,L)
call trajectory(r,N,time,counter)
call in_velocity(vel,N,Temp)
call rdf(r,N,L,0,nhis,density,delg,ngr,g)

!Cut-off calculation
cut=L*0.5d0
time=0.d0

!Initial forces,energies,pressure
call forces_LJ(L,N,r,cut,force,press,upot)
call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)

!Main Molecular Dynamics loop
do while (time.lt.Maxtime)
  time=time+dt
  counter=counter+1
  call verlet_velocity(N,cut,press,r,vel,force,dt,upot,L)
  call kinetic_en(vel,N,kin)
  call temperatura(kin,N,Temp)
  call trajectory(r,N,time,counter)
  call units_print(time,upot,kin,press,L,dt,sigma,eps,density,Temp,mass)
  call rdf(r,N,L,1,nhis,density,delg,ngr,g)
  call momentum(time,vel,N)
end do

call rdf(r,N,L,2,nhis,density,delg,ngr,g)

end program main
    
