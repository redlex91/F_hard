module utilities

  use prec_def
  use constants

  implicit none
  save
  real( pr ) :: sq_sigma, exp
  public :: sq_sigma, exp

contains

  subroutine init_cond( x, y, vx, vy )
    !Purpose : this sbrt initialises the system

    use prec_def
    use constants

    implicit none
    real( pr ), dimension( 1 : NOD ), intent( OUT ) ::  x, y, &
         vx, vy

    ! local variables
    integer :: i ! indexes for loops
    integer, parameter :: n = nint( sqrt( real( NOD, KIND = pr ) ) )
    real( pr ), parameter :: eps = 1.0_pr / real( n, KIND = pr )
    real( pr ) :: randnum, vxc, vyc, K
    integer, dimension( : ), allocatable :: seed
    integer, dimension( 1:8 ) :: dt_seed
    integer :: n_seed
    ! end local variables

    ! setting seed for speen generator
    call random_seed( size = n_seed )
    allocate( seed( 1 : n_seed ) )
    call random_seed( get = seed )
    call date_and_time( values = dt_seed )

    seed( n_seed ) = dt_seed( 8 ); seed( 1 ) = dt_seed( 8 ) * dt_seed( 7 ) * dt_seed( 6 )

    call random_seed( put = seed )
    deallocate( seed )
    ! done setting seed

    ! assigning initial coordinates
    do i = 1, NOD
       x( i ) = eps / 2._pr + mod( ( i - 1 ), n ) * eps 
       y( i ) = eps  / 2._pr + ( ( i - 1 ) / n ) * eps  
    end do

    ! assigning random velocities
    do i =1, NOD
       call random_number( randnum )
       vx( i ) = 2._pr * randnum - 1._pr
    end do
    do i =1 , NOD
       call random_number( randnum )
       vy( i ) = 2._pr * randnum - 1._pr
    end do

    ! computing the cm velocitiy
    vxc = 0._pr; vyc = 0._pr

    do i = 1, NOD
       vxc = vxc + vx( i ) / NOD
       vyc = vyc + vy( i ) / NOD
    end do
    ! subtracting the cm velocity
    do i = 1, NOD
       vx( i ) = vx( i ) -vxc
       vy( i ) = vy( i ) - vyc
    end do
    ! verifying cm velocity = 0
    vxc = 0._pr; vyc = 0._pr
    do i = 1, NOD
       vxc = vxc + vx( i )
       vyc = vyc + vy( i )
    end do
    ! PRINT *, vxc, vyc
    ! verified

    ! computing the kinetic energy...
    K = 0._pr
    !DO i = 1, NOD
    K = K + ( 1. / 2 ) * ( vx( i )**2 + vy( i )**2 )
    !END DO

    K = ( dot_product( vx, vx ) + dot_product( vy, vy ) ) / 2.0_pr 
    ! ... re-defining the velocities so that K = 1: this sets the scale for energy

    do i= 1 , NOD
       vx( i ) = vx( i ) / sqrt( K )
       vy( i ) = vy( i ) / sqrt( K )
    end do

  end subroutine init_cond ! end subroutine for defining the IC

  function sq_distance( x1, y1, x2, y2 )

    use prec_def

    implicit none
    real( pr ), intent( IN ) :: x1, x2, y1, y2
    real( pr ) :: sq_distance
    ! local var
    real( pr ) :: a1, a2, b1, b2
    a1 = x1; a2 = x2; b1 = y1; b2 = y2

    if( abs( a1 - a2 ) > 0.5_pr ) a2 = a2 - sign( 1.0_pr, a2 - a1 )
    if( abs( b1 - b2 ) > 0.5_pr ) b2 = b2 - sign( 1.0_pr, b2 - b1 )

    sq_distance = (a2 - a1 )**2._pr + ( b2 - b1 )**2._pr 

    return
  end function sq_distance

  function collision_time( x1, y1, vx1, vy1, x2, y2, vx2, vy2 )
    !Purpose: this sbrt  computes the collision time between two disks

    use prec_def
    use constants

    implicit none
    real( pr ), intent( IN ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2
    real( pr ) :: collision_time ! collision time
    ! local variables
    real( pr ), dimension( 1 : 9 ) :: vtimes ! this vector contains the collision times with the images
    real( pr ), dimension( 1 : 2 ) :: r12, v12
    real( pr ) :: r_dot_v, sqnr, sqnv, discrim ! dot product of r12 and v12, square module od r12, square module of v12
    integer :: p, q, k

    k = 1
    do p = -1 , 1
       do q = -1 , 1
          r12( 1 ) = x1 - ( x2 - p * 1._pr ); r12( 2 ) = y1 - ( y2 - q * 1._pr )
          v12( 1 ) = vx1 - vx2; v12( 2 ) = vy1 - vy2

          r_dot_v = dot_product( r12, v12 ) ! refer to formula 1.26
          sqnr = dot_product( r12, r12 )
          sqnv = dot_product( v12, v12 )

          if( sqnr - sq_sigma >= 0._pr ) then
             discrim = r_dot_v**(2._pr) - sqnv * ( sqnr - sq_sigma )
          else
             discrim = r_dot_v**2._pr
          end if

          ! IF(  sqnr - sigma**2._pr < 0._pr ) THEN
          !    PRINT *, "Error in collision time! Overlap: ", sqnr - sigma**2._pr
          !    CALL ABORT
          ! END IF

          if( r_dot_v < 0._pr .and. discrim >= 0._pr ) then ! there exists real solutions
             vtimes( k ) = ( - r_dot_v - sqrt( discrim ) ) / sqnv
             if( abs( vtimes( k ) ) < epsilon( 1._pr ) ) vtimes( k ) = 0._tpr 
          else
             vtimes( k ) = huge( 1.0_pr ) ! no collision: set tij to the maximum possible value
          end if
          k = k + 1
       end do
    end do

    ! PRINT *, "vtimes( k ): "
    ! DO k = 1, 9
    !    PRINT *, vtimes( k )
    ! END DO

    do k = 1, 9
       if( vtimes( k ) < 0. ) then
          print *, "ERROR in collision time! TIME:", vtimes( k )
          call ABORT
       end if
    end do

    collision_time = vtimes( 1 ) ! vtimes( 1 )
    do k = 1 , 9
       if ( vtimes( k ) < collision_time ) collision_time = vtimes( k )
    end do

    return
    !WRITE (*,*) "tempo di coll. ( compute_ct ): ", ctime, " dot p.:", r_dot_v

  end function collision_time


  ! SUBROUTINE ct_copies( x1, y1, vx1, vy1, x2, y2, vx2, vy2, mtime ) ! compute the collision time among a disks and the nine copies of a second

  !   USE prec_def
  !   USE constants

  !   IMPLICIT NONE
  !   REAL( pr ), INTENT( IN ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2
  !   REAL( pr ), INTENT( INOUT ) :: mtime ! minimum collision time among one disk and the nine copies of a second
  !   ! local var
  !   INTEGER :: p, q, k
  !   REAL( pr ), DIMENSION( 1:9 ) :: vtimes

  !   !WRITE (*,*) "COPIES"
  !   k=1
  !   DO p = -1,1
  !      DO q = -1,1
  !         !		WRITE (*,*), "vettore prima", vtimes(k)
  !         CALL compute_ct( x1, y1, vx1, vy1, x2 - p * 1._pr, y2 - q * 1._pr, vx2, vy2, vtimes( k ) )
  !         !		WRITE (*,*), "vettore dopo", vtimes(k)
  !         k = k + 1
  !      END DO
  !   END DO

  !   !DO k=1,9
  !   !	WRITE (*,*), "vettore dopo dopo", vtimes(k)
  !   !END DO

  !   CALL min_find( vtimes, mtime )

  !   !WRITE (*,*), "Min delle 9 copie: ", mtime
  ! END SUBROUTINE ct_copies


  ! SUBROUTINE min_find( vtimes, min_time )
  !   !Purpose: this sbrt finds the minimum of a vector

  !   USE prec_def
  !   USE constants

  !   IMPLICIT NONE
  !   REAL( pr ), DIMENSION( 1 : 9 ), INTENT( IN ) :: vtimes
  !   REAL( pr ), INTENT( INOUT ) :: min_time
  !   ! local variables
  !   INTEGER :: k

  !   min_time = HUGE( 1._pr ) ! vtimes( 1 )
  !   DO k = 1 , 9
  !      IF ( vtimes( k ) < min_time ) min_time = vtimes( k )
  !   END DO

  !   !WRITE (*,*), "Min vettore", min_time
  ! END SUBROUTINE min_find

  subroutine check_overlap( rx, ry, p1, p2 )
    use prec_def
    use constants

    implicit none
    real( pr ), intent( IN ) :: rx( 1 : NOD ), ry( 1 : NOD )
    integer, intent( IN ) :: p1, p2
    integer :: i, j

    do i = 1, NOD
       do j = i +  1, NOD
          if( ( ( i /= p1 .and. j /= p2 ) .and. ( j /= p1 .and. i /= p2 ) ) &
               .and. sq_distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sq_sigma < - 10 * epsilon( 1._pr ) ) then
             print *, 'The particles ', i, ' and  ', j, ' overlap!'
             print *, '1: ', rx( i ), ry( i ), '2: ', rx( j ), ry( j ), 'sigma: ', sqrt( sq_sigma )
             print *, 'The particles which collide are: ', p1, p2
             print *, sq_distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sq_sigma, ' smallest variation: ', epsilon( 1._pr )
             call ABORT
          end if
       end do
    end do
  end subroutine check_overlap

  subroutine check_time( table )
    use prec_def
    use constants

    implicit none
    real( tpr ), dimension( 1 : NOD, 1 : NOD ), intent( INOUT ) :: table
    integer :: i, j, err = 0, num = 0 ! num counts the elements in the table, err counts the huge elements in the table

    do i = 1, NOD
       do j = i + 1, NOD
          num = num + 1
          if( table( i, j ) < - epsilon( 1._pr ) ) then
             print *, 'There is a negative time: ', table( i, j )
             call ABORT
          elseif( abs( table( i, j ) - huge( 1._pr ) ) < epsilon( 1._pr ) ) then 
             err = err + 1
          end if
       end do
    end do
    if( err == num ) then
       print *, 'No collision'
       call abort
    end if

  end subroutine check_time

end module utilities
