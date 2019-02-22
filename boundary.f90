        SUBROUTINE BOUNDARY_CONDITIONS(r,n_part,L)
               IMPLICIT NONE
               REAL(8), INTENT(IN)  :: L
               REAL(8)              :: r(n_part,3)

               r = r-nint(r/L)*L
        END SUBROUTINE BOUNDARY_CONDITIONS

