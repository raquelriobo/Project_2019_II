subroutine boundary_conditions(r, n_part, L)
implicit none
integer :: i
integer, intent(in) :: n_part           ! Number of particles 
real(8), intent(in) :: L                ! Cell longitude
real(8)             :: r(n_part,3)      ! Coordinates matrix

! ### Computation of Boundary Conditions ### !
r =  r-nint(r/L)*L                            
end subroutine boundary_conditions
