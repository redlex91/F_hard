! Modules

!==============================================================================
!         START OF PROGRAM === START OF PROGRAM === START OF PROGRAM 
!==============================================================================

program hard_disks

  use prec_def
  use constants
  use utilities

  ! START OF PREAMBLE
  implicit none
  real( pr ), dimension( 1:NOD ) :: rx, ry, & ! positions
       vx, vy ! velocities

  real( tpr ), dimension( 1:NOD, 1:NOD ) :: time_tab ! collision times' table
  real( tpr ) :: t_ij ! collision time between the first colliding particles
  real( pr ) :: kin_en ! kinetic energy

  integer :: i, j, k, iter

  integer :: p1, p2 ! positions of colliding disks
  real( pr ) :: d
  real( tpr ) :: time = 0 ! time of the physical system: it increases by steps of t_ij
  real( pr) :: start_t, end_t

  ! eta variable
  real( pr ) :: ETA
  real(pr) :: ave_press = 0, un_press = 0 ! average pressure, uncertainty
  integer :: evar, meas_press = 0 ! percentage of eta, number of measures of the pressure during an integration period

  ! variables for updating speeds
  real(pr), dimension( 2 ) :: r1, r2, v1, v2, vers, v12 ! versor r12 = ( r2 - r1 ) / norm
  real(pr) :: norm ! norm( r2 - r1 )

  ! variables for computing pressure
  integer, parameter :: n_int = 300! number of step of integration + 2
  real( tpr ) :: t_int = 0
  real( pr ) :: sum_delta_v = 0, press = 0 ! press = pressure / ( NOD*T )
  real( pr ), dimension( 2 ) :: delta_v 

  ! variables for computing delta_r square
  integer :: flag = 0 ! it determines the first measure and the start of meas_time
  real( tpr ) :: meas_time = 0, start_time, step_time, diff_time
  real( pr ) :: xwrite, ywrite, msqd = 0 ! mean square displacement 
  integer :: sep,  max_sep = 0, sep_count ! separation between two measures of position, maximal separation, number of deltar for each separation
  integer :: eof = 0 ! it determines the end of a file
  real( pr ) :: stuff
  real( pr ), dimension( : , : ), allocatable :: posit_tab

  ! variables for computing the  mean free path and the mean free time
  real( pr ) :: mfp = 0, mft = 0
  integer :: count_coll = 0

  call cpu_time( start_t )
  open( unit = 4, file = "etavsp.dat", status = "replace", access = "sequential", position = "rewind" )   

  evar = 20
  eta_loop: DO
     
     select case ( evar )
     case ( 20 : 69 )
        evar = evar + 10
     case( 70 : 119 )
        evar = evar + 1
     case( 120 : 200 )
        evar = evar + 10
     case default
        exit
     end select
     
     ! FILES
     open( unit = 1, file = "energy.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of energies
     !open( unit = 3, file = "speed.dat", status = "unknown", access = "sequential", position = "append" )
     open( unit = 2, file = "pressure.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of pressures
     open( unit = 7, file = "deltarsq.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of deltarsq
     open( unit = 8, file = "times.dat", status = "replace", access = "sequential", position = "rewind" )
     open( unit = 13, file = "position.dat", status = "replace", access = "sequential", position = 'append' )

     !read( *, * ) evar
     !evar = 300.0
     count_coll = 0
     ave_press = 0
     un_press = 0
     meas_press = 0
     flag = 0

     ETA = ETA_M - ( evar  / 1000.0_pr ) * ETA_M
     print *, ETA

     sq_sigma = ( 4.0_pr * ETA ) / ( PI * NOD ) 
     forall( i = 1 : NOD ) ! initialising the table of collision times with all HUGEs
        forall( j = 1 : NOD )
           time_tab( i, j ) = huge( 1._pr )
        end forall
     end forall

     ! INITIALIZATION OF THE SYSTEM
     call init_cond( rx, ry, vx, vy ) 
     call check_overlap( rx, ry, 1, 2 )

     ! DO i = 1, NOD ! PRINTing the IC
     ! 	PRINT *, '(', rx( i ), ',', ry( i ), '),', '(', vx( i ), ',', vy( i ), ')'
     ! END DO

     time = 0
     ! EVOLUTION OF THE SYSTEM
     step: do iter = 1 , MAX_ITER ! one-step history

        time_update: if( iter /= 1 ) then
           do i = 1, NOD
              do j = i + 1, NOD
                 if( ( i == p1 .or. j == p2 .or. i == p2 .or. j == p1) ) then
                    ! we need to compute the collision time directly only for the collisions
                    ! involvin the particles which collided in the previous step
                    time_tab( i, j ) = collision_time( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ) )
                 end if
              end do
           end do
        else
           ! the table of collision times is generated for the first time
           do i= 1, NOD
              do j = i + 1 , NOD
                 time_tab( i, j ) = collision_time( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ) )
                 !PRINT *, i , j , time_tab( i, j )   
              end do
           end do
        end if time_update
        call check_time( time_tab )

        p1 = 1; p2 = 2
        t_ij = time_tab( 1, 2 ) ! computing the minimum time of the table
        compute_collision_time: do i = 1 , NOD
           do j = i + 1 , NOD
              if( time_tab( i, j ) < t_ij ) then
                 t_ij = time_tab( i, j ) 
                 p1 = i
                 p2 = j
              end if
           end do
        end do compute_collision_time

        count_coll = count_coll + 1

        ! print *, "Collision time: ", t_ij, "disks:", p1, p2

        ! BLOCK FOR MEASURES
        measures: if( iter >= N_TERM * NOD   ) then ! after termalisation
           first_meas: if( flag /=  0 ) then ! successive measures

              sum_delta_v = sqrt( dot_product( delta_v, delta_v ) ) + sum_delta_v  ! contribute to the virial (coming from the previous step)

              diff_time = meas_time + step_time - (time + t_ij )

              end_of_integ: if( diff_time < 0._pr ) then ! end of integration period
                 meas_time = meas_time + step_time
                 do i = 1, NOD
                    ! the system evolves by t_ij+diff_time < t_ij
                    rx( i ) = rx( i ) + ( t_ij + diff_time ) * vx( i )
                    ry( i ) = ry( i ) + ( t_ij + diff_time ) * vy( i )
                    rx( i ) = rx( i ) - floor( rx( i ) ) * 1._pr
                    ry( i ) = ry( i ) - floor( ry( i ) ) * 1._pr
                    ! storing the positions at regular intervals
                    ! for the mean square disp.
                    write( unit = 13, fmt = '(f9.6/f9.6)' ) rx( i ), ry( i )
                 end do

                 ! subtract elapsed time from the table of collision times
                 forall( i = 1 : NOD ) 
                    forall( j = i + 1 : NOD )
                       time_tab( i, j ) = time_tab( i, j ) - ( t_ij + diff_time )
                    end forall
                 end forall
                 call check_time( time_tab )
                 t_int = t_int + ( t_ij + diff_time )
                 time = time + (t_ij + diff_time )

                 ! storing the pressure
                 !  print *, t_int, step_time
                 kin_en = ( dot_product( vx, vx ) + dot_product( vy, vy ) ) / 2.0_pr
                 press = 1._pr + ( sqrt( sq_sigma ) / ( 3._pr * kin_en ) ) * ( sum_delta_v / t_int)
                 print *, t_int, step_time, press
                 write( unit = 2, fmt = '(f20.15)' ) press ! WRITE pressure data on file
                 ave_press = ave_press + press
                 meas_press = meas_press + 1 

                 ! refresh variables for new integration
                 sum_delta_v = 0
                 t_int = 0
              else ! no increment of meas_time during this step
                 diff_time = - t_ij
              end if end_of_integ

           else ! the integration starts for the first time
              print *, 'Starting measures.'
              meas_time = time
              start_time = time
              step_time = ceiling( 100 * 50._tpr * time / ( 2._tpr * count_coll / real( NOD, tpr ) ) )/100._tpr
              !WRITE( unit = 13, fmt = '(F8.6)' ) ( meas_time - start_time ) / step_time 
              do i = 1, NOD
                 write( unit = 13, fmt = '(F8.6/F8.6)' ) rx( i ), ry( i )
              end do
              !WRITE( unit = 13 , fmt = '()' )

              sum_delta_v = 0._pr
              t_int = 0._tpr
              diff_time = - t_ij
              flag = 1
           end if first_meas

           ! computing contributions to the mean free path
           do i = 1, NOD   
              mfp = mfp + sqrt( vx( i )**2._pr + vy( i )**2._pr ) * t_ij / NOD
           end do

           t_int = t_int + ( - diff_time )

        else ! before termalisation
           diff_time = - t_ij
        end if measures

        ! evolution of the system - new disks' positions
        forall( k = 1 : NOD ) ! assigns the new positions after the evolution
           rx( k ) = rx( k ) + vx( k ) * ( - diff_time )
           ry( k ) = ry( k ) + vy( k ) * ( - diff_time )
           rx( k ) = rx( k ) - floor( rx( k ) ) * 1._pr
           ry( k ) = ry( k ) - floor( ry( k ) ) * 1._pr
        end forall
        call check_overlap( rx, ry, p1, p2 )

        ! measuring total time
        time = time + ( - diff_time )

        !print *, time
        write( unit = 8, fmt = 200 ) t_ij, time
200     format( f20.8, f20.8 )

        ! EVOLUTION of the system - new disks' SPEEDS

        ! index of the vector: i=1 x-component, i=2 y-component
        r1( 1 ) = rx( p1 ); r1( 2 ) = ry( p1 )
        r2( 1 ) = rx( p2 ); r2( 2 ) = ry( p2 )
        v1( 1 ) = vx( p1 ); v1( 2 ) = vy( p1 )
        v2( 1 ) = vx( p2 ); v2( 2 ) = vy( p2 )

        ! chooses the minimum distance between two points in the box
        if( abs( r2(1) - r1(1) ) > 0.5_pr ) then
           r2(1) = r2(1) - sign( 1.0_pr, r2(1) - r1(1) )
        end if
        if( abs( r2(2) - r1(2) ) > 0.5_pr) then
           r2(2) = r2(2) - sign( 1.0_pr, r2(2) - r1(2) )
        end if

        norm = sqrt( dot_product( r2 - r1 , r2 - r1 ) )
        vers = ( r2 - r1 ) / norm 
        v12 = v1 - v2
        ! contribute to the virial of the colliding particles:
        ! delta_v = v1'-v1 - ( v2'-v2 ) = - dot_product( v12, vers ) * vers
        delta_v = - 2._pr * dot_product( v12 , vers ) * vers ! 

        v1 = v1 - dot_product( v12, vers ) * vers
        v2 = v2 + dot_product( v12, vers ) * vers

        ! assignment of the new speeds
        vx( p1 ) = v1( 1 ) ! first disk
        vy( p1 ) = v1( 2 ) 
        vx( p2 ) = v2( 1 ) ! second disk
        vy( p2 ) = v2( 2 ) 

        ! subtract elapsed time from the table of collision times
        forall( i = 1 : NOD )
           forall( j = i + 1 : NOD )
              time_tab( i, j ) = time_tab( i, j ) - (- diff_time ) 
           end forall
        end forall
        call check_time( time_tab )

        !write (*,400) ( iter*100.0_pr )/MAX_ITER
        !400     FORMAT("running...", f7.1, " %")

        !        PRINT *, time
        ! PRINT *, ( dot_PRODUCT( vx, vx ) + dot_PRODUCT( vy, vy ) ) / 2.0_pr
        write( unit = 1, fmt = '( f9.3 )' ) ( dot_product( vx, vx ) + dot_product( vy, vy ) ) / 2.0_pr ! writing energies on file

        ! ! PRESSURE COMPUTATION
        ! kin_en = ( dot_PRODUCT( vx, vx ) + dot_PRODUCT( vy, vy ) ) / 2.0_pr

        ! IF( iter >= N_TERM * NOD ) THEN ! the system is in equilibrium
        !    IF( MOD( iter, n_int ) /= 0 ) THEN 
        !       t_int =  t_int + t_ij
        !       sum_delta_v = SQRT( dot_PRODUCT( delta_v , delta_v ) ) + sum_delta_v
        !    ELSE
        !       IF( flag /= 0 ) THEN ! we discard the first value
        !          press = 1._pr + ( SQRT( sq_sigma ) / ( 3._pr * kin_en ) ) * ( sum_delta_v / t_int )
        !          WRITE( unit = 2, fmt = '(F20.15)' ) press ! WRITE pressure data on file
        !          ave_press = ave_press + press
        !          meas_press = meas_press + 1
        !       END IF
        !       flag = 1
        !       ! refresh variables for new integration
        !       t_int = 0
        !       sum_delta_v = 0
        !    END IF
        !    !PRINT *, iter, press, t_int
        ! END IF

        ! print speed at a time
        !  IF( iter == nspeed ) THEN
        !            DO i = 1, NOD
        !               WRITE( unit = 3, fmt = 700 ) vx(i)
        ! 700           FORMAT( f20.15 )
        !            END DO
        !         END IF


     end do step

     !=============================================================================
     !  END ITERATIONS === END ITERATIONS === END ITERATIONS === END ITERATIONS 
     !=============================================================================

     endfile 13

     ! computation of the mean square displacement
     rewind( unit =  13 )
     eof = 0;  k = 0
     do
        do i = 1, 2*NOD
           read( unit = 13, fmt = '(F8.6)',  IOSTAT = eof )  stuff
        end do
        if( eof /= 0 ) exit     
        k = k + 1
     end do
     !     PRINT *, k

     rewind 13
     max_sep = k
     print *, max_sep
     allocate( posit_tab( 1 : max_sep, 1 : 2*NOD ) )
     do k = 1, max_sep
        do i = 1, 2*NOD
           read( unit = 13, fmt = '(F8.6)' ) posit_tab( k, i )
        end do
     end do

     ! REWIND 13
     ! DO k = 1, max_sep
     !    DO i = 1, 2*NOD
     !       WRITE( unit = 6, fmt = '(F8.6)', ADVANCE = 'NO' )  posit_tab( k, i )
     !    END DO
     !    PRINT *, ' '
     ! END DO

     rewind 13
     do sep = 0, max_sep - 1
        msqd = 0; sep_count = 0
        do i = 1, max_sep
           do j = i, max_sep
              if( j - i == sep ) then
                 do k = 1, NOD
                    msqd = msqd + sq_distance( posit_tab( i, k ), posit_tab( i, 2*k ), posit_tab( j, k ), posit_tab( j, 2*k ) )
                 end do
                 sep_count = sep_count + 1
              end if
           end do
        end do
        msqd = msqd / real( NOD * sep_count )
        write( unit = 7, fmt = '(i5,F10.6 )' ) sep, msqd 
        ! PRINT *, sep_count
     end do

     deallocate( posit_tab )
     ! mean free path
     mfp = mfp / count_coll
     ! mean free time
     mft = time / ( ( 2._pr * count_coll ) / NOD )

     ! FILES
     close( unit = 1 )
     close( unit = 3 )
     close( unit = 7 )
     close( unit = 8 )
     close( unit = 13 )

     ! computing the average pressure as function of eta
     if( meas_press >= 2 ) then
        ave_press = ave_press / meas_press
        rewind 2
        eof = 0; k = 0; un_press = 0
        do
           read( unit = 2, fmt = '(F20.15)', IOSTAT = eof ) stuff
           un_press = un_press + ( stuff - ave_press )**2._pr
           if( eof /= 0 ) exit
           k = k + 1
        end do
        print *, k, ' = ', meas_press
        ! computing the uncertainty
        un_press = sqrt( un_press / ( meas_press * (meas_press - 1 ) ) )
        write( unit = 4, fmt = '( f10.3, f20.15, F20.15 )' ) eta, ave_press, un_press
     end if
     close( unit = 2, status = 'delete' )

  end do eta_loop

  close( unit = 4 )
  call cpu_time( end_t )
  write (*,900) end_t - start_t
900 format("Execution lasted: ", f10.2, " seconds.")

end program hard_disks

!==============================================================================
!            END OF PROGRAM === END OF PROGRAM === END OF PROGRAM
!==============================================================================

! the following is trash to be deleted!!!






   !SUBROUTINE evolve( rx, ry, vx, vy, ctime )
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
   !	rx( k ) = rx( k ) - FLOOR( rx( k ) )
   !	ry( k ) = ry( k ) - FLOOR( ry( k ) )
   !END DO

   !END SUBROUTINE evolve


   !SUBROUTINE up_speed( rx, ry, vx, vy, p1, p2 ) ! updates the speeds of the disks after the evolution
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
   !IF( ABS( r2(1) - r1(1) ) > 0.5 ) THEN
   !	r2(1) = r2(1) - SIGN( 1.0_prec, r2(1) - r1(1) )
   !END IF

   !IF( ABS( r2(2) - r1(2) ) > 0.5 ) THEN
   !	r2(2) = r2(2) - SIGN( 1.0_prec, r2(2) - r1(2) )
   !END IF


   ! WRITE (*,*), "Overlap?", distance( r1(1), r1(2), r2(1), r2(2) )-sigma

   !norm = sqrt( dot_PRODUCT( r2 - r1, r2 - r1 ) )
   !vers = ( r2 - r1 ) / norm 
   !v12 = v1 - v2

   !v1 = v1 - dot_PRODUCT( v12, vers ) * vers
   !v2 = v2 + dot_PRODUCT( v12, vers ) * vers

   ! assignment of the new speeds
   !vx( p1 ) = v1( 1 ) ! first disk
   !vy( p1 ) = v1( 2 ) 
   !vx( p2 ) = v2( 1 ) ! second disk
   !vy( p2 ) = v2( 2 ) 

   !end SUBROUTINE up_speed

   !     IF( iter < N_TERM * NOD ) THEN
   !        meas_time = time
   !     ELSEIF( iter == N_TERM * NOD ) THEN ! set instant t = 0 for the measure of the msqd
   !        meas_time = time + step_time
   !        FORALL( i = 1 : NOD )
   !           xd1( i ) = rx( i )
   !           yd1( i ) = ry( i )
   !           delta_r_sing( i ) = 0
   !        END FORALL
   ! !PRINT *, "PRIMA", meas_time
   !     END IF

   !      deltarsq: IF( meas_time - time <= t_ij .AND. iter >= N_TERM * NOD ) THEN
   ! 	delta1 = meas_time - time
   !  !PRINT *, iter, meas_time, delta1, t_ij
   !  ! evolution of the system by delta1
   ! 	DO k = 1 , NOD
   !            rx( k ) = rx( k ) + vx( k ) * delta1
   !            ry( k ) = ry( k ) + vy( k ) * delta1
   !            ! projection inside the box
   !            rx( k ) = rx( k ) - FLOOR( rx( k ) )
   !            ry( k ) = ry( k ) - FLOOR( ry( k ) )
   ! 	END DO
   !  ! measure of dr
   ! 	DO k = 1 , NOD
   !            delta_r_sing( k ) = distance( xd1( k ), yd1( k ), rx( k ), ry( k ) )**2._prec
   ! 	END DO
   ! 	delta_r = 0._prec
   ! 	DO i = 1 , NOD
   !            delta_r = delta_r + delta_r_sing( i )
   ! 	END DO
   ! 	WRITE( unit = 7, fmt = 100 ) delta_r/NOD
   ! 100     FORMAT( f20.15 )
   ! 	! translation of the table of collision times
   ! 	DO i = 1 , NOD
   !            DO j = i + 1 , NOD
   !               time_tab( i, j ) = time_tab( i, j ) - delta1
   !            END DO
   ! 	END DO
   !  ! refresh time of measure
   ! 	meas_time = meas_time + step_time
   !      ELSE 
   !         delta1 = 0
   !      END IF deltarsq

   !      delta2 = t_ij - delta1
