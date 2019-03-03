    subroutine boundary_conditions(r, n_part, L)
        implicit none
        integer :: i
	integer, intent(in) :: n_part 
        real(8), intent(in) :: L
        real(8)             :: r(n_part,3)
                
        r =  r-nint(r/L)*L
    end subroutine boundary_conditions
