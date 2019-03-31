subroutine temperatura(e_kin,n_part,temp)
implicit none
real(8), intent(in)    :: e_kin         ! Kinetic energy
real(8), intent(out)   :: temp          ! Temperature
integer, intent(in)    :: n_part        ! Number of particles

!Temperature calculation with equipartition theorem            
temp = (2.d0*e_kin)/(3.d0*n_part)

end subroutine temperatura
