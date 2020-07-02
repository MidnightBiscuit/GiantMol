:: This is a script batch for Windows
@echo off
Setlocal EnableDelayedExpansion
:: flags1    ='-xhost -O3 -align -save -qopenmp -mkl'
:: flags2='-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -save -qopenmp -par-num-threads=1  -mkl -warn unused'
set winflags1=/QxHost -O3 -align /Qsave /Qopenmp /Qmkl
set filename1=md_Langeving01
set filename_md=RF_relax_iter

echo =================
echo S
echo     T
echo         A
echo             R
echo                 T
echo =================

:: Udc
FOR %%e IN (9) DO (		
	set udc=000%%e
	set dir_name0=DC!udc:~-2!
	
	:: Vrf
	for %%g in (11) do (	
		set vrf=000%%g
		set dir_name1=RF!vrf:~-2!
		mkdir !dir_name0!_!dir_name1!
		
		:: try
		for /L %%t in (0,1,12) do (
			set T=000%%t
			set dir_name2=Try!T:~-2!
			
			set full_address=!dir_name0!_!dir_name1!\!dir_name2!
			echo !full_address!
			
			mkdir !full_address!
			echo =================
			echo   !full_address!
			echo =================

			:: Load from Langevin to md
			ifort -Dd0 -DnE00 -DHCI0 -DVrf%%g -DUdc%%e -DSimu_type3 %winflags1% %filename_md%.f90 -o !dir_name0!_!dir_name1!\!dir_name2!\a.out
			cd !dir_name0!_!dir_name1!\!dir_name2!
			a.out -ia 0 -ib 0
			cd ..\..
			
			:: Inject GMol
			ifort -Dd0 -DnE07 -DHCI1 -DVrf%%g -DUdc%%e -DSimu_type4 %winflags1% %filename_md%.f90 -o !dir_name0!_!dir_name1!\!dir_name2!\a.out
			cd !dir_name0!_!dir_name1!\!dir_name2!
			a.out -ia 0 -ib 0
			cd ..\..
			
			:: Post-injection dynamics
			ifort -Dd0 -DnE07 -DHCI0 -DVrf%%g -DUdc%%e -DSimu_type2 %winflags1% %filename_md%.f90 -o !dir_name0!_!dir_name1!\!dir_name2!\a.out
			cd !dir_name0!_!dir_name1!\!dir_name2!
			a.out -ia 0 -ib 0
			cd ..\..
		
	echo ======EGMol loop !udc:~-2!======
	echo N
	echo     E
	echo         X
	echo             T
	echo                 E
	echo                     G
	echo                         M
	echo ===========================
		)
	)
)