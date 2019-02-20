        SUBROUTINE BOUNDARY_CONDITIONS(n_part,r,L,dist_ij)
                IMPLICIT NONE
                
                INTEGER, INTENT(IN)  :: n_part
                REAL(8), INTENT(IN)  :: L, r(n_part,3)
                REAL(8), INTENT(OUT) :: dist_ij(3) 
                INTEGER              :: i,j,k


                do i=1,n_part
                  do j=1,n_part
                    if(j .ne. i) then
                        dist_ij(1) = r(i,1)-r(j,1)
                        dist_ij(2) = r(i,2)-r(j,2)
                        dist_ij(3) = r(i,3)-r(j,3)
                        do k=1,3
                            if(dist_ij(k) > L/2.) then
                                dist_ij(k) = dist_ij(k) - L
                            else if(dist_ij(k) <= -L/2.) then
                                dist_ij(k) = dist_ij(k) + L
                            end if
                        end do
                    end if 
                  end do
                end do

        END SUBROUTINE BOUNDARY_CONDITIONS

