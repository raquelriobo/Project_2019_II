! ### Binning ### !
program estadistica
implicit none
integer                             :: i, j, k
integer                             :: dim_binn, mbin, nbin
integer, parameter                  :: val_eliminat=0        ! Number of first deleted values
integer, parameter                  :: num_bin=11               ! Bin number
real(8), dimension(:), allocatable  :: epot, ekin, etot, p, temp
real(8)                             :: time, emitjana, varian

open(unit=25, file="Results.txt") 

! ###Subroutine that counts lines ### !
call count_lines(i)

allocate(epot(i))
allocate(ekin(i))
allocate(etot(i))
allocate(p(i))
allocate(temp(i))

! ### Read Results file and save values in different vectors ### !
do j=1, i
 read(25,*) time, epot(j), ekin(j), etot(j), p(j), temp(j)
end do
close(25)

! ### Subroutine Binning to different variables ###!

! ### New bin dimension without the first values ### !
dim_binn=i-val_eliminat

! ### Potential Energy ### !
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, epot(val_eliminat+1:i), varian, emitjana, 'Epotbinning.txt', 'Epotmit.txt')
end do

! ### Kinetic Energy ### !
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, ekin(val_eliminat+1:i), varian, emitjana, 'Ekinbinning.txt', 'Ekinmit.txt')
end do

! ### Total Energy ### !
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, etot(val_eliminat+1:i), varian, emitjana, 'Etotbinning.txt', 'Etotmit.txt')
end do

! ### Pressure ### !
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, p(val_eliminat+1:i), varian, emitjana, 'Presbinning.txt', 'Presmit.txt')
end do

! ### Temperature ### !
do k= 1, num_bin
 mbin=2**(k-1)
 nbin=int(dim_binn/mbin)
 call binning(dim_binn, mbin, nbin, temp(val_eliminat+1:i), varian, emitjana, 'Tempbinning.txt', 'Tempmit.txt')
end do

contains

! ### This subroutine count lines from a file and returns this value ### !
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

! ### Binning subroutine, prints in two different output files ### !
subroutine binning(dim_bin,mbin,nbin,variable,varian,emitjana,Filename,Namefile)
implicit none
character(15), intent(in) :: Filename !Filename for binning results
character(11), intent(in) :: Namefile !Filename for average value
integer:: dim_bin,mbin,nbin  
real(8):: promig(nbin), variable(dim_bin),var,varian,emitjana
integer:: i, j

open(unit=26, file=Filename)
open(unit=27, file=Namefile)

promig=0
do i=0,nbin-1
   do j=1,mbin
      promig(i+1)=promig(i+1)+variable(i*mbin+j)
   enddo
   promig(i+1)=promig(i+1)/float(mbin)
enddo

! ### Average Energy ### !
emitjana=0.0
do i=1,nbin
   emitjana=emitjana+promig(i)
end do
emitjana=emitjana/float(nbin)

! ### Variance ### !
var=0.0
do i=1,nbin
   var = var+(promig(i)-emitjana)**2
end do
var = var/(float(nbin)*(nbin-1))

varian=sqrt(var)

! ### Print results ###!
write(26,*) mbin, varian
write(27,*) emitjana
end subroutine

end program
