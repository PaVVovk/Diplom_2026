all: concentrations

# Исполняемый файл
concentrations: variables.o functions.o bolsig_data.o concentration.o electron_continuity.o ion_continuity.o Ambipolar.o BESSJ0.o Iterations.o Poisson.o Show_params.o
	gfortran variables.o functions.o bolsig_data.o concentration.o electron_continuity.o ion_continuity.o Ambipolar.o BESSJ0.o Iterations.o Poisson.o Show_params.o -o concentration.exe

# Объектные файлы
variables.o: variables.f90
	gfortran -Wall -c variables.f90

functions.o: functions.f90 variables.o bolsig_data.o
	gfortran -Wall -c functions.f90

bolsig_data.o: bolsig_data.f90 variables.o
	gfortran -Wall -c bolsig_data.f90

electron_continuity.o: electron_continuity.f90 variables.o functions.o
	gfortran -Wall -c electron_continuity.f90

ion_continuity.o: ion_continuity.f90 variables.o functions.o
	gfortran -Wall -c ion_continuity.f90

Ambipolar.o: Ambipolar.f90 variables.o
	gfortran -Wall -c Ambipolar.f90

BESSJ0.o: BESSJ0.f90
	gfortran -Wall -c BESSJ0.f90

Iterations.o: Iterations.f90 variables.o bolsig_data.o
	gfortran -Wall -c Iterations.f90

Poisson.o: Poisson.f90 variables.o 
	gfortran -Wall -c Poisson.f90

Show_params.o: Show_params.f90 variables.o functions.o
	gfortran -Wall -c Show_params.f90

concentration.o: concentration.f90 Show_params.o Iterations.o BESSJ0.o Ambipolar.o variables.o functions.o bolsig_data.o
	gfortran -Wall -c concentration.f90

clean:
	rm -f *.o *.exe

