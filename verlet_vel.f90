
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !          Each integration cycle:                 !
            !                                                  !
            !     CURRENT POS, VEL -> CURRENT FORCES           !
            !                                                  !
            !      CURRENT FORCES -> NEW POSITIONS :           !
            !                                                  !
            !    r(t+dt) = r(t) + v(t)*dt + F(t)/(2m)*dt^2     !
            !                                                  !
            !         NEW POSITIONS -> NEW FORCES              !
            !                                                  !
            !         NEW FORCES -> NEW VELOCITIES :           !
            !                                                  !
            !      v(t+dt) = v(t)+(F(t)+f(t+dt))/(2m)*dt       !
            !                                                  !  
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      SUBROUTINE VERLET_VELOCITY(n_part,cut_off,press,r,v,F,dt,E_pot,L)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n_part
            REAL(8), INTENT(IN) :: cut_off, dt
            REAL(8)             :: r(n_part,3), r_new(n_part,3) 
            REAL(8)             :: v(n_part,3), v_new(n_part,3)
            REAL(8)             :: F(n_part,3), F_new(n_part,3)
            REAL(8)             :: L,press, E_Pot
            INTEGER             :: i
              
              
            ! posicions i forces inicials
!            CALL forces_LJ_Press(L,n_part,r,cut_off,F,press,E_pot)

            ! posicions noves 
            do i=1,n_part
               r_new(i,:) = r(i,:) + v(i,:)*dt + &
               0.5*F(i,:)*dt**2
            end do 

            ! forces noves 
            CALL forces_LJ_Press(L,n_part,r_new,cut_off,F_new,press,E_pot)

            ! velocitats noves
            do i=1,n_part
                v_new(i,:) = v(i,:) + &
                0.5*(F(i,:)+F_new(i,:))*dt
            end do
                
            ! reinicialització 
            r = r_new ;    v = v_new ;    F = F_new ;
            r_new = 0.;    v_new = 0.;    F_new = 0.;
                   
	    call boundary_conditions(r,n_part,L) 
      END SUBROUTINE VERLET_VELOCITY
