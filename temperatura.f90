        SUBROUTINE TEMPERATURA(e_kin,n_part,temp)
            IMPLICIT NONE
            REAL, INTENT(IN)    :: e_kin
            REAL, INTENT(OUT)   :: temp
            INTEGER, INTENT(IN) :: n_part
            REAL                :: kb
            
            kb = 1.   !unitats reduides
            temp = (2.*e_kin)/(3.*n_part*kb)

        END SUBROUTINE TEMPERATURE
