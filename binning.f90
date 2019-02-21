integer,parameter:: mcs=1D7, res=1000, L=100    !Analisis de uno de los ficheros de salida de mc4f.90
real(8):: Epro, Mpro, Evar, Mvar,E(mcs),M(mcs)  !Parametros y vectores 
integer:: ii,aa,mbin,nbin,sd                     
sd=mcs-res                                      !dimension de la data para el binning (eliminacion de los 10^3 primeros)                          
open(8,file="temp3.0.txt")                      !Fichero de datos a leer
open(9,file="results3.0.txt")                   !Fichero  donde se guardan los resultados
do ii=1,mcs                                     !Bucle de lectura de los datos del fichero deseado
   read(8,*) a,E(ii),M(ii)                      
   M(ii)=abs(M(ii))                             !Pasamos la magnetizacion a valor absoluto
enddo
close(8)                                        !Cerramos el fichero de lectura
do ii=1,22                                      !Bucle para distintas m en el binning
   mbin=2**(ii-1)                               !Valor de m
   nbin=int(sd/mbin)                            !Valor de la longitud del vector del bining
   CALL binning(sd,mbin,nbin,E(res+1:mcs),Evar,Epro)   !Binning energia eliminando los 1000 primeros valores 
   CALL binning(sd,mbin,nbin,M(res+1:mcs),Mvar,Mpro)   !Binning magnetizacion eliminando los 1000 primeros valores
   write(9,"(I8,F15.10,F15.10,F15.10,F15.10)") mbin,Epro/(L**2),&  !Escritura de los resultados divididos entre el 
   Evar/(L**2),Mpro/(L**2),Mvar/(L**2)                             !numero de spines ya que eran energias totales
enddo
close(9)                                        !Cerramos el archivo de los resultados


contains


SUBROUTINE binning(n,m,mb,v,rf,Vm)    !Subrutina para hacer el bining dado un vector de dimension n,
integer:: n,m,mb                      !el valor de la m y el vector. Devuelve rf (eerror estandar) y
real(8):: pro(mb), v(n),var,rf,Vm     !Vm (valor promedio)
integer:: ii, jj
pro=0
do ii=0,mb-1
   do jj=1,m
      pro(ii+1)=pro(ii+1)+v(ii*m+jj)
   enddo
   pro(ii+1)=pro(ii+1)/float(m)
enddo
!!valor medio
Vm=0
do ii=1,mb
   Vm=Vm+pro(ii)
enddo
Vm=Vm/float(mb)
!! variança
var=0
do ii=1,mb
   var=var+(pro(ii)-Vm)**2
enddo
var=var/(float(mb)*(mb-1))
!!desviacio
rf=sqrt(var)
end subroutine


end program
