program estadistica
implicit none
integer :: i, j, k, dim_binn, mbin, nbin
real*8, dimension(:), allocatable :: epot, ekin, etot, p, temp
real*8 :: time, emitjana, varian
character(15) :: Epotbinning, Ekinbinning, Etotbinning
integer, parameter :: val_eliminat=100! número de dades que s'eliminen del principi
integer, parameter :: num_bin=10 !!poden ser més si cal
open(unit=25, file="Results.txt") 

call count_lines(i)

allocate(epot(i))
allocate(ekin(i))
allocate(etot(i))
allocate(p(i))
allocate(temp(i))

do j=1, i
 read(25,*) time, epot(j), ekin(j), etot(j), p(j), temp(j)
end do
close(25)
!Aplica la subrutina de binning i escriu els resultats en els arxius
!corresponents
dim_binn=i-val_eliminat
!Energia potencial
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, epot(val_eliminat+1:i), varian, emitjana, 'Epotbinning.txt')
end do
!Energia cinètica
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, ekin(val_eliminat+1:i), varian, emitjana, 'Ekinbinning.txt')
end do
!Energia total
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, etot(val_eliminat+1:i), varian, emitjana, 'Etotbinning.txt')
end do
!Pressió
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, epot(val_eliminat+1:i), varian, emitjana, 'Presbinning.txt')
end do
!Temperatura
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, epot(val_eliminat+1:i), varian, emitjana, 'Tempbinning.txt')
end do

contains

subroutine count_lines(num_lines)
implicit none
integer, intent(out) :: num_lines
open(unit=25, file="Results.txt")

num_lines=0
do
 read(25, *, END=52)
 num_lines=num_lines+1
end do
52 rewind(25)

end subroutine

subroutine binning(dim_bin,mbin,nbin,variable,varian,emitjana,Filename)
implicit none
character(15), intent(in) :: Filename
integer:: dim_bin,mbin,nbin  
real(8):: promig(nbin), variable(dim_bin),var,varian,emitjana
integer:: i, j

open(unit=26, file=Filename)
promig=0
do i=0,nbin-1
   do j=1,mbin
      promig(i+1)=promig(i+1)+variable(i*mbin+j)
   enddo
   promig(i+1)=promig(i+1)/float(mbin)
enddo
! Calcula l'Energia mitjana del sistema
emitjana=0.0
do i=1,nbin
   emitjana=emitjana+promig(i)
end do
emitjana=emitjana/float(nbin)

! Calcula la variança del sistema
var=0.0
do i=1,nbin
   var = var+(promig(i)-emitjana)**2
end do
var = var/(float(nbin)*(nbin-1))

!Desviació estàndard d'aquestes dades
varian=sqrt(var)

!imprimeix els resultats a l'arxiu corresponent
write(26,*) mbin, emitjana, varian
end subroutine

end program
