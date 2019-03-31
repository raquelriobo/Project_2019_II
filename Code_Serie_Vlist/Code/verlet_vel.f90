! ### INTEGRATION : VERLET VELOCITY #### !
subroutine verlet_velocity(n_part,cut_off,press,r,v,F,dt,E_pot,L,nlist,list,contF,timeforces)
implicit none
integer, intent(in) :: n_part                       ! Number of particles
integer :: contF
real(8), intent(in) :: cut_off                      
real(8), intent(in) :: dt                           ! Time steps
real(8)             :: r(n_part,3), r_new(n_part,3) ! Coordinates matrix 
real(8)             :: v(n_part,3), v_new(n_part,3) ! Velocity matrix
real(8)             :: F(n_part,3), F_new(n_part,3) ! Forces matrix
real(8)             :: L                            ! Cell longitude
real(8)             :: press                        ! Pressure
real(8)             :: E_Pot                        ! Potential energy
integer             :: i
!real(8) :: startF(150000), finishF(150000)
real(8) :: timeforces(150000)
integer :: nlist(n_part),list(n_part,n_part-1)              
              
contF=contF+1
! ### Compute new positions with verlet velocity algorithm ###! 
do i=1,n_part
    r_new(i,:) = r(i,:) + v(i,:)*dt + &
    0.5*F(i,:)*dt**2
end do 

! ### Compute new forces ###!
CALL forces_vlist(L,n_part,r_new,cut_off,F_new,press,E_pot,nlist,list,contF,timeforces)

!print*, contF

! ### Compute new velocity ###! 
do i=1,n_part
    v_new(i,:) = v(i,:) + &
    0.5*(F(i,:)+F_new(i,:))*dt
end do
                
! ### Update #### ! 
r = r_new ;    v = v_new ;    F = F_new ;
r_new = 0.;    v_new = 0.;    F_new = 0.;

call boundary_conditions(r,n_part,L) 

end subroutine verlet_velocity
