        SUBROUTINE TEMPERATURA(e_kin,n_part,temp)
            IMPLICIT NONE
            REAL(8), INTENT(IN)    :: e_kin
            REAL(8), INTENT(OUT)   :: temp
            INTEGER, INTENT(IN)    :: n_part
            REAL(8)                :: kb
            
            kb = 1.d0   !unitats reduides
            temp = (2.d0*e_kin)/(3.d0*n_part*kb)

        END SUBROUTINE TEMPERATURA
