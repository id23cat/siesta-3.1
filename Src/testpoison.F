	program testpoison
C		poison( CELL, N1, N2, N3, Mesh, RHO, U, V, STRESS, 
C     &                   NSM )
C**************** INPUT **********************************************
C REAL*8  CELL(3,3)     : UNIT CELL VECTORS
C INTEGER N1,N2,N3      : NUMBER OF MESH DIVISIONS IN EACH CELL VECTOR
C INTEGER Mesh(3)       : Number of global mesh divisions
C REAL*4  RHO(N1,N2,N3) : DENSITIY AT MESH POINTS
C **************** OUTPUT *********************************************
C REAL*8  U             : ELECTROSTATIC ENERGY (IN RY)
C REAL*4  V(N1,N2,N3)   : ELECTROSTATIC POTENTIAL (IN RY)
C                         V AND RHO MAY BE THE SAME PHYSICAL ARRAY
C REAL*8  STRESS(3,3) : Electrostatic-energy contribution to stress
C                       tensor (in Ry/Bohr**3) assuming constant density
C                       (not charge), i.e. r->r' => rho'(r') = rho(r)
C                       For plane-wave and grid (finite difference)
C                       basis sets, density rescaling gives an extra
C                       term (not included) equal to -2*U/cell_volume
C                       for the diagonal elements of stress. For other
C                       basis sets, the extra term is, in general:
C                       IntegralOf( V * d_rho/d_strain ) / cell_volume
C INTEGER NSM           : Number of sub-mesh points per mesh point
C                       : along each axis
!	use	m_dhscf, only : dhscf_init
	use precision,	only : dp, grid_p
	!use m_timer, only: timer_init, timer_start, timer_stop, timer_report
	use mesh,	only : nsm
!	use cudafor	!only : cudaEvent_t, cudaEventCreate, cudaEventDestroy, & 
			!cudaEventRecord, cudaEventSynchronize, cudaEventElapsedTime
	
	!include 'curand.h'
	!external nsm

	
	!curandGenerator_t	gen
	
	implicit none

	integer,parameter	:: N=500
	integer,parameter	:: LOOPS=1
	real(dp)		:: CELL(3,3), STRESS_CPU(3,3), STRESS_GPU(3,3), UC, UG
	integer,parameter	:: N1=N, N2=N, N3=N
	
	integer			:: mesh(3)
	real(grid_p)		:: RHO(N1,N2,N3), VC(N1*N2*N3), VG(N1*N2*N3)
	integer			:: i,j,k, istat

!	type(cudaEvent)		:: startEv, stopEv;
	integer 		:: c0, c1, c2, c3, clock_rate
	real			:: cgpu, chost
	
	
	
	mesh(1) = N1
	mesh(2) = N2
	mesh(3) = N3
	!nsm = fdf_integer( 'MeshSubDivisions', 2 )
	nsm = 2
	do i=1,3
	do j=1,3
		CELL(i,j) = (i-1)*3 + j;
	end do
	end do
	
	do i=1,N1
	do j=1,N2
	do k=1,N3
		RHO(i,j,k) = 11.99+i+j+k !N1*N2*N3/((i-1)*N1+(j-1)*N2+1)
	end do
	end do
	end do
	
!	istat = cudaEventCreate(startEv)
!	istat = cudaEventCreate(stopEv)
	
!	call timer_init()
!	call timer_start('ACC Poison')
!	istat = cudaEventRecord(startEv,0)

	CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
	call  SYSTEM_CLOCK( COUNT=c1 )
	
	call poison(CELL, N1, N2, N3, Mesh, RHO, UG, VG, STRESS_GPU, nsm)
	
	call  SYSTEM_CLOCK( COUNT=c2 )
!	istat = cudaEventRecord(stopEv, 0)
!	istat = cudaEventElapsedTime(elapsedTime, startEv, stopEv)
!	write(*,*) 'ACC - elpsed time:\t\t', elapsedTime
!	call timer_stop('ACC Poison')
	!call timer_report(printNow=.true.)
	
!	call timer_start('Original Poison')
!	istat = cudaEventRecord(startEv,0)
	
	call poison_orig(CELL, N1, N2, N3, Mesh, RHO, UC, VC, STRESS_GPU, nsm)
	
	call  SYSTEM_CLOCK( COUNT=c3 )
	
	cgpu=REAL((c2-c1)/clock_rate)
	cgpu=REAL((c3-c2)/clock_rate)
	
	print *, cgpu, ' microseconds on GPU'
	print *, chost, ' microseconds on host'
	
!	istat = cudaEventRecord(stopEv, 0)
!	istat = cudaEventElapsedTime(elapsedTime, startEv, stopEv)
!	write(*,*) 'Original - elpsed time:\t', elapsedTime
!	call timer_stop('Original Poison')
	!call timer_report(printNow=.true.)
	
!	do i=1,3
!	do j=1,3
!		if(STRESS_CPU(i,j).ne.STRESS_GPU(i,j))
!			print *, ' unequal in ', i,j,k
!	end do
!	end do

!	write(*,*) STRESS_CPU
!	write(*,*) STRESS_GPU
!	write(*,*) 'U'
!	write(*,*) UC
!	write(*,*) UG
!	write(*,*) 'V'
!	write(*,*) VC
!	write(*,*) VG

!	write(*,*) RHO
	
	end program testpoison
