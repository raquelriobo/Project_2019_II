#Makefile

##backup : Copia de seguridad comprimida
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90


