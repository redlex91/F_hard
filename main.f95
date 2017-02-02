! modules

MODULE prec_def
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER:: sp = KIND( 1.0 ), &
       dp = SELECTED_REAL_KIND( 2 * PRECISION( 1._sp ) ), &
       qp = SELECTED_REAL_KIND( 2 * PRECISION( 1._dp ) ), &
       prec = qp ! precision in use
ENDMODULE prec_def

MODULE constants
  USE prec_def

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: NOD = 100, N_TERM = 13 ! NOD = number of disks
  INTEGER, PARAMETER :: MAX_ITER = 150000, nspeed = 30 ! max number of iterations

  REAL( prec ), PARAMETER :: PI = 3.14159265358979323846264338327950288419716939937510_prec, ETA_M = PI/4
ENDMODULE constants

MODULE utilities

  USE prec_def
  USE constants
  
  IMPLICIT NONE
  SAVE
  REAL( prec ) :: sigma
  PUBLIC :: sigma

CONTAINS
  
  SUBROUTINE init_cond( x, y, vx, vy )
    !Purpose : this sbrt initialises the system
    
    USE prec_def
    USE constants

    IMPLICIT NONE
    REAL( prec ), DIMENSION( 1 : NOD ), INTENT( OUT ) ::  x, y, &
         vx, vy

    ! local variables
    INTEGER :: i ! indexes for loops
    INTEGER, PARAMETER :: n = NINT( SQRT( REAL( NOD ) ) )
    REAL( prec ), PARAMETER :: eps = 1.0 / n
    REAL( prec ) :: randnum, vxc, vyc, K
    INTEGER, DIMENSION( : ), ALLOCATABLE :: seed
    INTEGER, DIMENSION( 1:8 ) :: dt_seed
    INTEGER :: n_seed
    ! end local variables

    ! setting seed for speen generator
    CALL random_seed( size = n_seed )
    ALLOCATE( seed( 1:n_seed ) )
    CALL random_seed( get = seed )
    CALL date_and_time( values = dt_seed )

    seed( n_seed ) = dt_seed( 8 ); seed( 1 ) = dt_seed( 8 ) * dt_seed( 7 ) * dt_seed( 6 )

    CALL random_seed( put = seed )
    DEALLOCATE( seed )
    ! done setting seed

    ! assigning initial coordinates
    DO i = 1 , NOD
       x( i ) = eps / 2 + MOD( ( i - 1 ), n ) * eps 
       y( i ) = eps  / 2 + ( ( i - 1 ) / n ) * eps  
    ENDDO

    ! assigning random velocities
    DO i =1 , NOD
       CALL random_number( randnum )
       vx( i ) = 2._prec * randnum - 1._prec
    ENDDO
    DO i =1 , NOD
       CALL random_number( randnum )
       vy( i ) = 2._prec * randnum - 1._prec
    ENDDO

    ! computing the cm velocitiy
    vxc = 0._prec; vyc = 0._prec

    DO i = 1, NOD
       vxc = vxc + vx( i ) / NOD
       vyc = vyc + vy( i ) / NOD
    ENDDO
    ! subtracting the cm velocity
    DO i = 1, NOD
       vx( i ) = vx( i ) -vxc
       vy( i ) = vy( i ) - vyc
    ENDDO
    ! verifying cm velocity = 0
    vxc = 0; vyc = 0
    DO i = 1, NOD
       vxc = vxc + vx( i )
       vyc = vyc + vy( i )
    ENDDO
    ! print *, vxc, vyc
    ! verified

    ! computing the kinetic energy...
    K = 0
    !DO i = 1, NOD
    K = K + ( 1. / 2 ) * ( vx( i )**2 + vy( i )**2 )
    !ENDDO

    K = ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) ) / 2.0_prec 
    ! ... re-defining the velocities so that K = 1: this sets the scale for energy

    DO i= 1 , NOD
       vx( i ) = vx( i ) / SQRT( K )
       vy( i ) = vy( i ) / SQRT( K )
    ENDDO

  ENDSUBROUTINE init_cond ! end subroutine for defining the IC

  FUNCTION distance( x1, y1, x2, y2 )

    USE prec_def

    IMPLICIT NONE
    REAL( prec ), INTENT( IN ) :: x1, x2, y1, y2
    REAL( prec ) :: distance
    ! local var
    REAL( prec ) :: a1, a2, b1, b2
    a1 = x1; a2 = x2; b1 = y1; b2 = y2

    IF( ABS( a1 - a2 ) > 0.5_prec ) a2 = a2 - SIGN( 1.0_prec, a2 - a1 )
    IF( ABS( b1 - b2 ) > 0.5_prec ) b2 = b2 - SIGN( 1.0_prec, b2 - b1 )

    distance = SQRT( (a2 - a1 )**2._prec + ( b2 - b1 )**2._prec )
   
    RETURN
  ENDFUNCTION distance

  FUNCTION collision_time( x1, y1, vx1, vy1, x2, y2, vx2, vy2 )
    !Purpose: this sbrt  computes the collision time between two disks

    USE prec_def
    USE constants

    IMPLICIT NONE
    REAL( prec ), INTENT( IN ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2
    REAL( prec ) :: collision_time ! collision time
    ! local variables
    REAL( prec ), DIMENSION( 1 : 9 ) :: vtimes ! this vector contains the collision times with the images
    REAL( prec ), DIMENSION( 1 : 2 ) :: r12, v12
    REAL( prec ) :: r_dot_v, sqnr, sqnv, discrim ! dot product of r12 and v12, square module od r12, square module of v12
    INTEGER :: p, q, k

    k = 1
    DO p = -1 , 1
       DO q = -1 , 1
          r12( 1 ) = x1 - ( x2 - p * 1._prec ); r12( 2 ) = y1 - ( y2 - q * 1._prec )
          v12( 1 ) = vx1 - vx2; v12( 2 ) = vy1 - vy2

          r_dot_v = DOT_PRODUCT( r12, v12 ) ! refer to formula 1.26
          sqnr = DOT_PRODUCT( r12, r12 )
          sqnv = DOT_PRODUCT( v12, v12 )
          discrim = r_dot_v**2._prec - sqnv * ( sqnr - sigma**2._prec )

          IF( r_dot_v < 0._prec .AND. discrim >= 0._prec ) THEN ! there exists real solutions
             vtimes( k ) = ( - r_dot_v - SQRT( discrim ) ) / sqnv 
          ELSE
             vtimes( k ) = HUGE( 1.0_prec ) ! no collision: set tij to the maximum possible value
          ENDIF
          k = k + 1
       END DO
    END DO

    collision_time = HUGE( 1._prec ) ! vtimes( 1 )
    DO k = 1 , 9
       IF ( vtimes( k ) < collision_time ) collision_time = vtimes( k )
    ENDDO


    !write (*,*) "tempo di coll. ( compute_ct ): ", ctime, " dot p.:", r_dot_v

  ENDFUNCTION collision_time
  
  ! SUBROUTINE ct_copies( x1, y1, vx1, vy1, x2, y2, vx2, vy2, mtime ) ! compute the collision time among a disks and the nine copies of a second

  !   USE prec_def
  !   USE constants

  !   IMPLICIT NONE
  !   REAL( prec ), INTENT( IN ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2
  !   REAL( prec ), INTENT( INOUT ) :: mtime ! minimum collision time among one disk and the nine copies of a second
  !   ! local var
  !   INTEGER :: p, q, k
  !   REAL( prec ), DIMENSION( 1:9 ) :: vtimes

  !   !write (*,*) "COPIES"
  !   k=1
  !   DO p = -1,1
  !      DO q = -1,1
  !         !		write (*,*), "vettore prima", vtimes(k)
  !         CALL compute_ct( x1, y1, vx1, vy1, x2 - p * 1._prec, y2 - q * 1._prec, vx2, vy2, vtimes( k ) )
  !         !		write (*,*), "vettore dopo", vtimes(k)
  !         k = k + 1
  !      ENDDO
  !   ENDDO

  !   !DO k=1,9
  !   !	write (*,*), "vettore dopo dopo", vtimes(k)
  !   !endDO

  !   CALL min_find( vtimes, mtime )

  !   !write (*,*), "Min delle 9 copie: ", mtime
  ! ENDSUBROUTINE ct_copies

  
  ! SUBROUTINE min_find( vtimes, min_time )
  !   !Purpose: this sbrt finds the minimum of a vector

  !   USE prec_def
  !   USE constants

  !   IMPLICIT NONE
  !   REAL( prec ), DIMENSION( 1 : 9 ), INTENT( IN ) :: vtimes
  !   REAL( prec ), INTENT( INOUT ) :: min_time
  !   ! local variables
  !   INTEGER :: k

  !   min_time = HUGE( 1._prec ) ! vtimes( 1 )
  !   DO k = 1 , 9
  !      IF ( vtimes( k ) < min_time ) min_time = vtimes( k )
  !   ENDDO

  !   !write (*,*), "Min vettore", min_time
  ! ENDSUBROUTINE min_find

ENDMODULE utilities












































PROGRAM hard_disks

  USE prec_def
  USE constants
  USE utilities
  
  ! DICHIARAZIONE DELLE VARIABILI
  IMPLICIT NONE
  REAL( prec ), DIMENSION( 1:NOD ) :: rx, ry, & ! positions
       vx, vy ! velocities

  REAL( prec ), DIMENSION( 1:NOD, 1:NOD ) :: time_tab ! collision times' table
  REAL( prec ) :: t_ij ! collision time between the first colliding particles
  REAL( prec ) :: kin_en ! kinetic energy

  INTEGER :: i, j, k, iter

  INTEGER :: p1, p2 ! positions of colliding disks
  REAL( prec ) :: d
  REAL( prec) :: start_t, end_t

  ! eta variable
  REAL( prec ) :: ETA
  REAL(prec) :: meanP = 0
  INTEGER :: evar, Pdata = 0

  ! variables for updating speeds
  REAL(prec), DIMENSION( 2 ) :: r1, r2, v1, v2, vers, v12 ! versor r12 = ( r2 - r1 ) / norm
  REAL(prec) :: norm ! norm( r2 - r1 )

  ! variables for computing pressure
  INTEGER, PARAMETER :: n_int = 30! number of step of integration + 2
  REAL( prec ) :: t_int = 0, sum_delta_v = 0, press = 0 ! press = pressure / ( NOD*T )
  REAL( prec ), DIMENSION( 2 ) :: delta_v 
  INTEGER :: flag = 0


  ! variables for computing delta_r square
  REAL( prec ), DIMENSION( 1:NOD ) :: xd1, yd1 !, yd1, yd2 ! positions of particle at start time and end time
  REAL( prec ), DIMENSION( 1:NOD ) :: delta_r_sing 
  REAL( prec ) :: delta_r
  REAL( prec ) :: meas_time = 0, step_time = 0.10, time = 0, delta1, delta2
  INTEGER :: count_dr, flag_dr


  ! FILES
  open( unit = 1, file = "energy.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of energies
  open( unit = 2, file = "pressure.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of pressures
  open( unit = 3, file = "speed.dat", status = "unknown", access = "sequential", position = "append" )
  open( unit = 4, file = "etavsp.dat", status = "unknown", access = "sequential", position = "append" )
  open( unit = 7, file = "deltarsq.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of deltarsq
  open( unit = 8, file = "times.dat", status = "replace", access = "sequential", position = "rewind" )

  call cpu_time( start_t )


  !read( *, * ) evar
  evar = 500.0
  ETA = ETA_M - ( evar  / 1000.0_prec ) * ETA_M
  !print *, ETA

  sigma = sqrt( ( 4.0_prec * ETA ) / ( PI * NOD ) )

  CALL init_cond( rx, ry, vx, vy ) ! initiliasing the system

  ! DO i = 1, NOD ! printing the IC
  ! 	print *, '(', rx( i ), ',', ry( i ), '),', '(', vx( i ), ',', vy( i ), ')'
  ! ENDDO

  FORALL( i = 1 : NOD ) ! initialising the table of collision times with all 0s
     FORALL( j = 1 : NOD )
        time_tab( i, j ) = 0._prec
     ENDFORALL
  ENDFORALL


  ! I write the delta_r squared for t = 0: delta_r = 0
  write( unit = 7, fmt = 101 ) 0.0
101 format( f20.5 )
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  ! 																ITERATIONS === ITERATIONS === ITERATIONS === ITERATIONS === ITERATIONS 
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================

  PRINT *, "sigma = ", sigma
  step: DO iter = 1 , MAX_ITER ! one-step history



     !write (*,*), ""
     !write (*,*), "PASSO", iter
     !write (*,*), ""
     ! calcolo energia
     !write (*,*), "energia: ",

     ! verifies that no overlap has happened
     !DO i=1,NOD 
     !	DO j=1, NOD
     !		if( distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma < -10**(-10) .and. ( i .ne. j ) ) then
     !			print *, "ERRORE, passo", iter
     !			call abort
     !		endif
     !	ENDDO
     !ENDDO


     time_refresh: IF( iter == 1 ) THEN ! first step
        DO i= 1, NOD
           DO j = i , NOD
              IF( i == j ) THEN
                 time_tab( i, j ) = 0
              ELSE
                 time_tab( i, j ) = collision_time( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ) )
                 !PRINT *, i , j , time_tab( i, j )
              ENDIF
           ENDDO
        ENDDO
     ELSE ! following steps
        DO i = 1 , NOD
           DO j = i + 1 , NOD
              IF( ( i == p1 .or. j == p2 .or. i == p2 .or. j == p1) ) THEN ! interaction with one of the two last-colliding disks
                 time_tab( i, j ) = collision_time( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ) )
              ELSE
                 ! other disks
                 time_tab( i, j ) = time_tab( i, j ) - t_ij ! - delta2
              ENDIF
           ENDDO
        ENDDO
     ENDIF time_refresh


     	! DO i=1, NOD
     	! 	DO j=i,NOD
     	! 		write (*, *) time_tab( i, j )
     	! 	endDO
     	! 	print *, ""
     	! endDO

     t_ij = HUGE( 1.0_prec ) ! computing the minimum time of the table
     compute_collision_time: DO i = 1 , NOD
        DO j = i , NOD
           IF( time_tab( i, j ) < t_ij .AND. i.NE.j ) THEN
              t_ij = time_tab( i, j ) 
              p1 = i
              p2 = j
           ENDIF
        ENDDO
     ENDDO compute_collision_time
     !print *, ""

     IF( t_ij < 0 ) THEN
	PRINT *, " tempo <0 "
	CALL abort
     ENDIF
     print *, "Collision time: ", t_ij, "disks:", p1, p2

 !     IF( iter < N_TERM * NOD ) THEN
 !        meas_time = time
 !     ELSEIF( iter == N_TERM * NOD ) THEN ! set instant t = 0 for the measure of the msqd
 !        meas_time = time + step_time
 !        FORALL( i = 1 : NOD )
 !           xd1( i ) = rx( i )
 !           yd1( i ) = ry( i )
 !           delta_r_sing( i ) = 0
 !        ENDFORALL
 ! !print *, "PRIMA", meas_time
 !     ENDIF

!      deltarsq: IF( meas_time - time <= t_ij .AND. iter >= N_TERM * NOD ) THEN
! 	delta1 = meas_time - time
!  !print *, iter, meas_time, delta1, t_ij
!  ! evolution of the system by delta1
! 	DO k = 1 , NOD
!            rx( k ) = rx( k ) + vx( k ) * delta1
!            ry( k ) = ry( k ) + vy( k ) * delta1
!            ! projection inside the box
!            rx( k ) = rx( k ) - FLOOR( rx( k ) )
!            ry( k ) = ry( k ) - FLOOR( ry( k ) )
! 	ENDDO
!  ! measure of dr
! 	DO k = 1 , NOD
!            delta_r_sing( k ) = distance( xd1( k ), yd1( k ), rx( k ), ry( k ) )**2._prec
! 	ENDDO
! 	delta_r = 0._prec
! 	DO i = 1 , NOD
!            delta_r = delta_r + delta_r_sing( i )
! 	ENDDO
! 	WRITE( unit = 7, fmt = 100 ) delta_r/NOD
! 100     FORMAT( f20.15 )
! 	! translation of the table of collision times
! 	DO i = 1 , NOD
!            DO j = i + 1 , NOD
!               time_tab( i, j ) = time_tab( i, j ) - delta1
!            ENDDO
! 	ENDDO
!  ! refresh time of measure
! 	meas_time = meas_time + step_time
!      ELSE 
!         delta1 = 0
!      ENDIF deltarsq

!      delta2 = t_ij - delta1

     ! evolution of the system - new disks' positions
     FORALL( k = 1 : NOD ) ! assigns the new positions after the evolution
	rx( k ) = rx( k ) + vx( k ) * t_ij!delta2
	ry( k ) = ry( k ) + vy( k ) * t_ij!delta2
 ! projection inside the box
	rx( k ) = rx( k ) - FLOOR( rx( k ) )
	ry( k ) = ry( k ) - FLOOR( ry( k ) )
     ENDFORALL


     ! measuring total time
     time = time + t_ij
     !print *, time
     WRITE( unit = 8, fmt = 200 ) t_ij, time
200  FORMAT( f20.15, f20.15 )


     ! EVOLUTION of the system - new disks' SPEEDS
     !call up_speed( rx, ry, vx, vy, p1, p2 )

     ! index of the vector: i=1 x-component, i=2 y-component
     r1( 1 ) = rx( p1 ); r1( 2 ) = ry( p1 )
     r2( 1 ) = rx( p2 ); r2( 2 ) = ry( p2 )
     v1( 1 ) = vx( p1 ); v1( 2 ) = vy( p1 )
     v2( 1 ) = vx( p2 ); v2( 2 ) = vy( p2 )

     ! chooses the minimum distance between two points in the box
     IF( ABS( r2(1) - r1(1) ) > 0.5_prec ) THEN
	r2(1) = r2(1) - SIGN( 1.0_prec, r2(1) - r1(1) )
     ENDIF

     IF( ABS( r2(2) - r1(2) ) > 0.5_prec) THEN
	r2(2) = r2(2) - SIGN( 1.0_prec, r2(2) - r1(2) )
     ENDIF

     ! write (*,*), "Overlap?", distance( r1(1), r1(2), r2(1), r2(2) )-sigma

     norm = SQRT( DOT_PRODUCT( r2 - r1 , r2 - r1 ) )
     vers = ( r2 - r1 ) / norm 
     v12 = v1 - v2

     delta_v = 2 * DOT_PRODUCT( v12 , vers ) * vers 

     v1 = v1 - DOT_PRODUCT( v12, vers ) * vers
     v2 = v2 + DOT_PRODUCT( v12, vers ) * vers

     ! assignment of the new speeds
     vx( p1 ) = v1( 1 ) ! first disk
     vy( p1 ) = v1( 2 ) 
     vx( p2 ) = v2( 1 ) ! second disk
     vy( p2 ) = v2( 2 ) 


     !print *, ""
     !DO i = 1, NOD ! printing the IC
     !	print *, '(', rx( i ), ',', ry( i ), '),', '(', vx( i ), ',', vy( i ), ')'
     !ENDDO

     !write (*,*), "energia: ",( dot_product( vx, vx ) + dot_product( vy, vy ) )/2.0_prec
     ! OVERLAP ROUTINE : verifies that no overlap has happened
     DO i = 1 , NOD 
        DO j = i + 1 , NOD
           IF( distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma  < -TINY( 1._prec )  &
                .AND.  i .NE. p1  .AND.  j .NE. p2  )  THEN
              PRINT *, "ERRORE, passo", iter, distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma, "disks: ", i, j
              CALL abort
           ENDIF
        ENDDO
     ENDDO



     WRITE (*,400) ( iter*100.0_prec )/MAX_ITER
400  FORMAT("running...", f7.1, " %")

     WRITE( unit = 1, fmt = 500 ) ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) ) / 2.0_prec ! writing energies on file
500  FORMAT( f9.3 )

   !   ! PRESSURE COMPUTATION
!      kin_en = ( DOT_PRODUCT( vx, vx ) + DOT_PRODUCT( vy, vy ) ) / 2.0_prec
!      IF( MOD( iter, n_int ) == 0 .AND. iter >= N_TERM * NOD ) THEN 
! 	IF( flag .NE. 0 ) THEN
!            press = 1._prec + ( ( sigma ) / ( 3._prec * kin_en ) ) * ( sum_delta_v / t_int )
!            WRITE( unit = 2, fmt = 600 ) press ! write pressure data on file
! 600        FORMAT( f20.15 )
!            meanP = meanP + press
!            Pdata = Pdata + 1
! 	ENDIF
!  !print *, iter, press, t_int
! 	t_int = 0
! 	sum_delta_v = 0
!      ELSEIF( ( MOD( iter, n_int ) .NE. 0 ) .AND. ( iter >= N_TERM * NOD ) ) THEN 
! 	flag = 1
! 	t_int =  t_int + t_ij
! 	sum_delta_v = SQRT( DOT_PRODUCT( delta_v , delta_v ) ) + sum_delta_v
!      ENDIF

     ! print speed at a time
     IF( iter == nspeed ) THEN
	DO i = 1, NOD
           WRITE( unit = 3, fmt = 700 ) vx(i)
700        FORMAT( f20.15 )
	ENDDO
     ENDIF

  ENDDO step

  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  ! 														END ITERATIONS === END ITERATIONS === END ITERATIONS === END ITERATIONS === END ITERATIONS 
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================
  !===========================================================================================================================================================================================

  !meanP = meanP / Pdata
  WRITE( unit = 4, fmt = 800 ) eta, meanP
800 FORMAT( f10.3, f20.15 )

  CALL cpu_time( end_t )
  WRITE (*,900) end_t - start_t
900 FORMAT("Execution lasted: ", f10.2, " seconds.")

  ! FILES
  CLOSE( unit = 1 )
  CLOSE( unit = 2 )
  CLOSE( unit = 3 )
  CLOSE( unit = 4 )
  CLOSE( unit = 7 )
  CLOSE( unit = 8 )

ENDPROGRAM hard_disks




!===========================================================================================================================================================================================
!================================================================================== SUBROUTINE DEFINITIONS =================================================================================
!===========================================================================================================================================================================================




!subroutine evolve( rx, ry, vx, vy, ctime )
!USE prec_def
!USE constants
!IMPLICIT NONE
!REAL( prec ), DIMENSION( 1:NOD ), intent( in ) :: vx, vy
!REAL( prec ), intent( in ) :: ctime
!REAL( prec ), DIMENSION( 1:NOD ), intent( inout ) :: rx, ry
! local var
!INTEGER k

!DO k = 1, NOD
!	rx( k ) = rx( k ) + vx( k ) * ctime
!	ry( k ) = ry( k ) + vy( k ) * ctime
!	! projection inside the box
!	rx( k ) = rx( k ) - floor( rx( k ) )
!	ry( k ) = ry( k ) - floor( ry( k ) )
!ENDDO

!end subroutine evolve


!subroutine up_speed( rx, ry, vx, vy, p1, p2 ) ! updates the speeds of the disks after the evolution
!USE prec_def
!USE constants
!IMPLICIT NONE
!REAL( prec ), DIMENSION( 1:NOD ), intent( in ) :: rx, ry
!REAL( prec ), DIMENSION( 1:NOD ), intent( inout ) :: vx, vy
!INTEGER, intent( in ) :: p1, p2
! local variables
!REAL(prec), DIMENSION( 2 ) :: r1, r2, v1, v2, vers, v12 ! versor r12 = ( r2 - r1 ) / norm
!REAL(prec) :: norm, d, distance ! norm( r2 - r1 )

! index of the vector: i=1 x-component, i=2 y-component
!r1( 1 ) = rx( p1 ); r1( 2 ) = ry( p1 )
!r2( 1 ) = rx( p2 ); r2( 2 ) = ry( p2 )
!v1( 1 ) = vx( p1 ); v1( 2 ) = vy( p1 )
!v2( 1 ) = vx( p2 ); v2( 2 ) = vy( p2 )

! chooses the minimum distance between two points in the box
!if( abs( r2(1) - r1(1) ) > 0.5 ) then
!	r2(1) = r2(1) - sign( 1.0_prec, r2(1) - r1(1) )
!endif

!if( abs( r2(2) - r1(2) ) > 0.5 ) then
!	r2(2) = r2(2) - sign( 1.0_prec, r2(2) - r1(2) )
!endif


! write (*,*), "Overlap?", distance( r1(1), r1(2), r2(1), r2(2) )-sigma

!norm = sqrt( dot_product( r2 - r1, r2 - r1 ) )
!vers = ( r2 - r1 ) / norm 
!v12 = v1 - v2

!v1 = v1 - dot_product( v12, vers ) * vers
!v2 = v2 + dot_product( v12, vers ) * vers

! assignment of the new speeds
!vx( p1 ) = v1( 1 ) ! first disk
!vy( p1 ) = v1( 2 ) 
!vx( p2 ) = v2( 1 ) ! second disk
!vy( p2 ) = v2( 2 ) 

!endsubroutine up_speed
