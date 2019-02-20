
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



        SUBROUTINE VERLET_VELOCITY(n_part,r_ini,v_ini,mass,eps,sigma)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: n_part
                REAL(8), INTENT(IN) :: r_ini(n_part,3)
                REAL(8), INTENT(IN) :: v_ini(n_part,3)
                REAL(8), INTENT(IN) :: eps,sigma
                REAL(8)             :: r(n_part,3), r_new(n_part,3) 
                REAL(8)             :: v(n_part,3), v_new(n_part,3)
                REAL(8)             :: F(n_part,3), F_new(n_part,3)
                REAL(8)             :: mass, Energia_cin, U_LJ
                REAL(8)             :: ek,dist_ij(3),dx,dy,dz,d2, &
                                       cut_off
                INTEGER             :: i,j,iter
                INTEGER, PARAMETER  :: n_iteracions = 10
                REAL, PARAMETER     :: dt = 0.1
                
                
                !Inicialització valors posició, velocitat, força
                r = r_ini
                v = v_ini
                ! F_inicial == call Forces(..)

                do iter = 1, n_iteracions

                    ! og.forces -> new positions
                    do i=1,n_part
                        r_new(i,:) = r(i,:) + v(i,:)*dt + &
                        F(i,:)/(2.*mass)*dt**2
                    end do 

                    ! new positions -> new forces

                    !CALL FORCES(...)

                    ! new forces -> new velocities
                    do i=1,n_part
                        v_new(i,:) = v(i,:) + &
                        (F(i,:)+F_new(i,:))/(2.*mass)*dt
                    end do
                

                    ! CALCUL ENERGIA CINETICA
                    do i=1,n_part
                        ek = 0.5*mass*DOT_PRODUCT(v_new(i,:),v_new(i,:))
                        Energia_cin = Energia_cin + ek
                    end do

                    ! CALCUL ENERGIA POTENCIAL
                    do i=1,n_part
                        dist_ij(1) = dx(r_new(i,1),r_new(j,1))
                        dist_ij(2) = dy(r_new(i,2),r_new(j,2))
                        dist_ij(3) = dz(r_new(i,3),r_new(j,3))
                        d2 = DOT_PRODUCT(dist_ij,dist_ij)

                        if(sqrt(d2)<=cut_off .and. sqrt(d2)>0.00001) &
                        then
                            U_LJ = U_LJ + 2.*eps* & 
                            (sigma**12/d2**6 - sigma**6/d2**3)
                        end if
                    end do

                    ! reinicialització per executar nova iteració
                    r = r_new
                    r_new = 0.
                    
                    v = v_new
                    v_new = 0.
                    
                    F = F_new
                    F_new = 0.

                end do

        END SUBROUTINE VERLET_VELOCITY
