! *** Total momentum *** !
SUBROUTINE momentum(time,vel,N,p)
implicit none
integer :: N,i
real(8) :: time,vel(N,3),p(3)
open(15,file="v_verlet_momnt.txt")
p=0d0

do i=1,N ! * Loop over particles
	p(1)=p(1)+vel(i,1) ! * px
	p(2)=p(2)+vel(i,2) ! * py
	p(3)=p(3)+vel(i,3) ! * pz
end do

WRITE(15,*) time,p(1),p(2),p(3)

END SUBROUTINE momentum

