! modules
module prec_def
implicit none
save
integer, parameter :: long = selected_real_kind( 18,307 )
end module prec_def

module constants
use prec_def
implicit none
save
integer, parameter :: NOD = 225, N_TERM = 7 ! NOD = number of disks
real(long), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_long, ETA_M = PI/4
end module constants

program hard_disks

use prec_def
use constants

! DICHIARAZIONE DELLE VARIABILI
implicit none
real(long), dimension( 1:NOD ) :: rx, ry, vx, vy
real( long ), dimension( 1:NOD, 1:NOD ) :: time_tab ! collision times' table
real( long ) :: t_ij ! collision time between the first colliding particles
real( long ) :: kin_en ! kinetic energy
integer :: i, j, k, iter
integer, parameter :: MAX_ITER = 1e5, nspeed = 30 ! max number of iterations
integer :: p1, p2 ! positions of colliding disks
real( long ) :: d, distance ! function00
real( long) :: start_t, end_t

! eta variable
real( long ) :: ETA
real(long) :: sigma, meanP = 0
integer :: evar, Pdata = 0

! variables for updating speeds
real(long), dimension( 2 ) :: r1, r2, v1, v2, vers, v12 ! versor r12 = ( r2 - r1 ) / norm
real(long) :: norm ! norm( r2 - r1 )

! variables for computing pressure
integer, parameter :: n_int = 30! number of step of integration + 2
real( long ) :: t_int = 0, sum_delta_v = 0, press = 0 ! press = pressure / ( NOD*T )
real( long ), dimension( 2 ) :: delta_v 
integer :: flag = 0


! variables for computing delta_r square
real( long ), dimension( 1:NOD ) :: xd1, xd2, yd1, yd2 ! positions of particle at start time and end time
real( long ), dimension( 1:NOD ) :: delta_r, delta_r_sing 
real( long ) :: start_time, end_time, step_time, chrono = 0, evo_time
integer :: count_dr, flag_dr

interface ! prototypes for variables check

subroutine init_cond( x, y, vx, vy ) ! assigns the initial conditions fo the system
use prec_def
use constants
implicit none
real( long ), dimension( 1 : NOD ), intent( out ) ::  x, y, vx, vy
end subroutine init_cond

subroutine compute_ct( x1, y1, vx1, vy1, x2, y2, vx2, vy2, time, sigma ) ! compute the collision time between two disks
use prec_def
use constants
implicit none
real( long ), intent( in ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2, sigma
real( long ), intent( inout ) :: time
endsubroutine compute_ct

subroutine min_find( vtimes, min_time ) ! finds the minimum of a vector
use prec_def
use constants
implicit none
real( long ), dimension( 1:9 ), intent( in ) :: vtimes
real( long ), intent( inout ) :: min_time
endsubroutine min_find

subroutine ct_copies( x1, y1, vx1, vy1, x2, y2, vx2, vy2, mtime, sigma ) ! compute the collision time between a disk 
! and one of the nine copies of another
use prec_def
use constants
implicit none
real( long ), intent( in ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2, sigma
real( long ), intent( inout ) :: mtime ! minimum collision time among one disk and the nine copies of a second
endsubroutine ct_copies

!subroutine evolve( rx, ry, vx, vy, ctime ) ! assigns new values to r after the evolution 
!use prec_def
!use constants
!implicit none
!real( long ), dimension( 1:NOD ), intent( in ) :: vx, vy
!real( long ), intent( in ) :: ctime
!real( long ), dimension( 1:NOD ), intent( inout ) :: rx, ry
!endsubroutine evolve

!subroutine up_speed( rx, ry, vx, vy, p1, p2 )
!use prec_def
!use constants
!implicit none
!real( long ), dimension( 1:NOD ), intent( in ) :: rx, ry
!real( long ), dimension( 1:NOD ), intent( inout ) :: vx, vy
!integer, intent( in ) :: p1, p2
!endsubroutine up_speed

end interface ! end prototypes for variable check

call cpu_time( start_t )


read( *, * ) evar

ETA = ETA_M - ( evar  / 1000.0_long ) * ETA_M
print *, ETA

sigma = sqrt( ( 4.0 * ETA ) / ( PI * NOD ) )

call init_cond( rx, ry, vx, vy ) ! initiliasing the system

!do i = 1, NOD ! printing the IC
!	print *, '(', rx( i ), ',', ry( i ), '),', '(', vx( i ), ',', vy( i ), ')'
!end do

do i = 1, NOD ! initialising the table of collision times with all 0s
	do j= 1,NOD
		time_tab( i, j ) = 0
	end do
end do

! FILES
open( unit = 1, file = "energy.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of energies
open( unit = 2, file = "pressure.dat", status = "replace", access = "sequential", position = "rewind" ) ! opening file of pressures
open( unit = 3, file = "speed.dat", status = "unknown", access = "sequential", position = "append" )
open( unit = 4, file = "etavsp.dat", status = "unknown", access = "sequential", position = "append" )
open( unit = 5, file = "deltarsq.dat", status = "replace", access = "sequential", position = "append" ) ! opening file of deltarsq
!===========================================================================================================================================================================================
!===========================================================================================================================================================================================
!===========================================================================================================================================================================================
! 																ITERATIONS === ITERATIONS === ITERATIONS === ITERATIONS === ITERATIONS 
!===========================================================================================================================================================================================
!===========================================================================================================================================================================================
!===========================================================================================================================================================================================

step: do iter = 1, MAX_ITER ! one-step history

	
	
!write (*,*), ""
!write (*,*), "PASSO", iter
!write (*,*), ""
! calcolo energia
!write (*,*), "energia: ",

! verifies that no overlap has happened
!do i=1,NOD 
!	do j=1, NOD
!		if( distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma < -10**(-10) .and. ( i .ne. j ) ) then
!			print *, "ERRORE, passo", iter
!			call abort
!		endif
!	end do
!end do


	time_ref: if( iter == 1 ) then ! first step
		do i= 1, NOD
			do j= i, NOD
				if( i == j ) then
					time_tab( i, j ) = 0
				else
					call ct_copies( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ), time_tab( i, j ), sigma )
				end if
			end do
		end do
	else ! following steps
		do i = 1, NOD
			do j = i + 1, NOD
				if( ( i == p1 .or. j == p2 .or. i == p2 .or. j == p1) ) then ! interaction with one of the two last-colliding disks
					call ct_copies( rx( i ), ry( i ), vx( i ), vy( i ), rx( j ), ry( j ), vx( j ), vy( j ), time_tab( i, j ), sigma )
				else ! other disks
					time_tab( i, j ) = time_tab( i, j ) - t_ij
				end if
			end do
		end do
	end if time_ref

!	do i=1, NOD
!		do j=i,NOD
!			write (*, *), time_tab( i, j )
!		enddo
!		print *, ""
!	enddo
	
	t_ij = huge( 1.0_long ) ! computing the minimum time of the table
	collision_time: do i = 1, NOD
		do j = i , NOD
			if ( time_tab( i, j ) < t_ij .and. i.ne.j ) then
				t_ij = time_tab( i, j ) 
				p1 = i
				p2 = j
			end if
		end do
	end do collision_time
!print *, ""

if( t_ij < 0 ) call abort
!print *, "Collision time: ", t_ij, "disks:", p1, p2

! measuring total time
chrono = chrono + t_ij
print *, chrono

! start of delta_r_sq

start_dr: if( start_time >= t_ij ) then

! I'm computing the time between the last collsion and start_time
	evo_time = t_ij - ( chrono - start_time )
	! I let the system evolve for evo_time
	do i = 1, NOD
		xd1( i ) = rx( i ) + vx( i ) * evo_time
		yd1( i ) = ry( i ) + vy( i ) * evo_time
		! I project the new positions if the disks go outside the box
			xd1( i ) = xd1( i ) - floor( xd1( i ) )
			yd1( i ) = yd1( i ) - floor( yd1( i ) )
		endif
	enddo
	! updating start_time
	start_time = start_time + step_time
endif start_dr
end_dr: if( end_time > t_ij ) then

	evo_time = t_ij - ( chrono - t_end )
	 ! I let the system evolve for evo_time
	do i = 1, NOD
		xd2( i ) = rx( i ) + vx( i ) * evo_time
		yd2( i ) = ry( i ) + vy( i ) * evo_time
		! I project the new positions if the disks go outside the box
			xd2( i ) = xd1( i ) - floor( xd1( i ) )
			yd2( i ) = yd1( i ) - floor( yd1( i ) )
		endif
	enddo
	! updating start_time
	end_time = end_time + step_time
	! I compute the delta_r for a single particle
	delta_r
	
endif end_dr
! end of delta_r_sq


! evolution of the system - new disks' positions
do k = 1, NOD ! assigns the new positions after the evolution
	rx( k ) = rx( k ) + vx( k ) * t_ij
	ry( k ) = ry( k ) + vy( k ) * t_ij
	! projection inside the box
	rx( k ) = rx( k ) - floor( rx( k ) )
	ry( k ) = ry( k ) - floor( ry( k ) )
end do

! EVOLUTION of the system - new disks' SPEEDS
!call up_speed( rx, ry, vx, vy, p1, p2 )

! index of the vector: i=1 x-component, i=2 y-component
r1( 1 ) = rx( p1 ); r1( 2 ) = ry( p1 )
r2( 1 ) = rx( p2 ); r2( 2 ) = ry( p2 )
v1( 1 ) = vx( p1 ); v1( 2 ) = vy( p1 )
v2( 1 ) = vx( p2 ); v2( 2 ) = vy( p2 )

! chooses the minimum distance between two points in the box
if( abs( r2(1) - r1(1) ) > 0.5 ) then
	r2(1) = r2(1) - sign( 1.0_long, r2(1) - r1(1) )
endif

if( abs( r2(2) - r1(2) ) > 0.5 ) then
	r2(2) = r2(2) - sign( 1.0_long, r2(2) - r1(2) )
endif

! write (*,*), "Overlap?", distance( r1(1), r1(2), r2(1), r2(2) )-sigma

norm = sqrt( dot_product( r2 - r1, r2 - r1 ) )
vers = ( r2 - r1 ) / norm 
v12 = v1 - v2

delta_v = 2 * dot_product( v12, vers ) * vers

v1 = v1 - dot_product( v12, vers ) * vers
v2 = v2 + dot_product( v12, vers ) * vers

! assignment of the new speeds
vx( p1 ) = v1( 1 ) ! first disk
vy( p1 ) = v1( 2 ) 
vx( p2 ) = v2( 1 ) ! second disk
vy( p2 ) = v2( 2 ) 


!print *, ""
!do i = 1, NOD ! printing the IC
!	print *, '(', rx( i ), ',', ry( i ), '),', '(', vx( i ), ',', vy( i ), ')'
!end do

!write (*,*), "energia: ",( dot_product( vx, vx ) + dot_product( vy, vy ) )/2.0_long
! OVERLAP ROUTINE : verifies that no overlap has happened
do i=1,NOD 
	do j= i + 1 , NOD
		if( ( distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma  < -1E-13 ) .and. ( i .ne. p1 ) .and. ( j .ne. p2 )) then
			print *, "ERRORE, passo", iter, distance( rx( i ), ry( i ), rx( j ), ry( j ) ) - sigma, "disks: ", i, j
			call abort
		endif
	end do
end do
write (*,100) ( iter*100.0_long )/MAX_ITER
	100 format ("running...", f7.1, " %")

write( unit = 1, fmt = 11 ) ( dot_product( vx, vx ) + dot_product( vy, vy ) )/2.0_long ! writing energies on file
	11 format( f20.18 )
	
! PRESSURE COMPUTATION
kin_en = ( dot_product( vx, vx ) + dot_product( vy, vy ) )/2.0_long
if( mod( iter, n_int ) == 0 .and. iter >= N_TERM * NOD ) then 
	if( flag .ne. 0 ) then
		press = 1 + ( ( sigma ) / ( 3 * kin_en ) ) * ( sum_delta_v / t_int )
		write( unit = 2, fmt = 12 ) press ! write pressure data on file
			12 format( f20.15 )
			meanP = meanP + press
			Pdata = Pdata + 1
	endif
	!print *, iter, press, t_int
	t_int = 0
	sum_delta_v = 0
else if( ( mod( iter, n_int ).ne.0 ) .and. ( iter >= N_TERM * NOD ) ) then 
	flag = 1
	t_int =  t_int + t_ij
	sum_delta_v = sqrt( dot_product( delta_v, delta_v ) ) + sum_delta_v
endif

! print speed at a time
if( iter == nspeed ) then
	do i = 1, NOD
		write( unit = 3, fmt = 103 ) vx(i)
			103 format( f20.15 )
	enddo
endif

end do step


meanP = meanP / Pdata
write( unit = 4, fmt = 33 ) eta, meanP
	33 format( f10.3, f20.15 )

call cpu_time( end_t )
write (*,110) end_t - start_t
	110 format ("Execution lasted: ", f10.2, " seconds.")
	
! FILES
close( unit = 1 )
close( unit = 2 )
close( unit = 3 )
close( unit = 4 )
close( unit = 5 )
endprogram hard_disks









subroutine init_cond( x, y, vx, vy )
use prec_def
use constants

implicit none
real( long ), dimension( 1 : NOD ), intent( out ) ::  x, y, vx, vy

! local variables
integer :: i ! indexes for cycles
integer, parameter :: n = nint( sqrt( real ( NOD ) ) )
real( long ), parameter :: eps = 1.0 / n
real( long ) :: randnum, vxc, vyc, K
integer, dimension( : ), allocatable:: seed
integer, dimension( 1:8 ) :: dt_seed
integer :: n_seed
! end local variables

! setting seed for speen generator
call random_seed( size = n_seed )
allocate( seed( 1:n_seed ) )
call random_seed( get = seed )
call date_and_time( values = dt_seed )
seed( n_seed ) = dt_seed( 8 ); seed( 1 ) = dt_seed( 8 ) * dt_seed( 7 ) * dt_seed( 6 )
call random_seed( put = seed )
deallocate( seed )
! done setting seed

! assigning initial coordinates
do i = 1, NOD
	x( i ) = eps / 2 + mod( ( i - 1 ), n ) * eps 
	y( i ) = eps  / 2 + ( ( i - 1 ) / n ) * eps  
end do

! assigning random velocities
do i =1, NOD
	call random_number( randnum )
	vx( i ) = 2 * randnum -1
end do
do i =1, NOD
	call random_number( randnum )
	vy( i ) = 2 * randnum -1
end do

! computing the cm velocitiy
vxc = 0; vyc = 0
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
vxc = 0; vyc = 0
do i = 1, NOD
	vxc = vxc + vx( i )
	vyc = vyc + vy( i )
end do
! print *, vxc, vyc
! verified

! computing the kinetic energy...
K = 0
!do i = 1, NOD
	K = K + ( 1. / 2 ) * ( vx( i )**2 + vy( i )**2 )
!end do

K = (dot_product( vx, vx ) + dot_product( vy, vy ) )/2.0_long 
! ... re-defining the velocities so that K = 1
do i= 1, NOD
	vx( i ) = vx( i ) / sqrt( K )
	vy( i ) = vy( i ) / sqrt( K )
end do

end subroutine init_cond ! end subroutine for defining the IC

subroutine compute_ct( x1, y1, vx1, vy1, x2, y2, vx2, vy2, time, sigma ) ! compute the collision time between two disks
use prec_def
use constants
implicit none
real( long ), intent( in ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2, sigma
real( long ), intent( inout ) :: time
! local variables
real( long ), dimension( 1:2 ) :: r12, v12
real( long ) :: dp, sqnr, sqnv, discrim

r12( 1 ) = x1 - x2; r12( 2 ) = y1 - y2
v12( 1 ) = vx1 - vx2; v12( 2 ) = vy1 - vy2

dp = dot_product( r12, v12 )
sqnr = dot_product( r12, r12 )
sqnv = dot_product( v12, v12 )
discrim = dp*dp - sqnv * ( sqnr - sigma * sigma )

if( dp <0 .and. discrim >= 0 ) then
	time = ( - dp - sqrt( discrim ) ) / sqnv 
	else
		time = huge( 1.0_long )
endif

! write (*,*), "tempo di coll.: ", time, " dot p.:", dp

endsubroutine compute_ct

subroutine min_find( vtimes, min_time ) ! finds the minimum of a vector
use prec_def
use constants
implicit none
real( long ), dimension( 1:9 ), intent( in ) :: vtimes
real( long ), intent( inout ) :: min_time
! local variables
integer :: k

min_time = vtimes( 1 )
do k = 1,9
	if ( vtimes( k ) < min_time ) min_time = vtimes( k )
end do

!write (*,*), "Min vettore", min_time
endsubroutine min_find

subroutine ct_copies( x1, y1, vx1, vy1, x2, y2, vx2, vy2, mtime, sigma ) ! compute the collision time among a disks and the nine copies of a second
use prec_def
use constants
implicit none
real( long ), intent( in ) :: x1, y1, vx1, vy1, x2, y2, vx2, vy2, sigma
real( long ), intent( inout ) :: mtime ! minimum collision time among one disk and the nine copies of a second
! local var
integer :: p, q, k
real( long ), dimension( 1:9 ) :: vtimes

!write (*,*) "COPIES"
	

k=1
do p = -1,1
	do q = -1,1
!		write (*,*), "vettore prima", vtimes(k)
		call compute_ct( x1, y1, vx1, vy1, x2 - p*1.0, y2 - q*1.0, vx2, vy2, vtimes( k ), sigma )
!		write (*,*), "vettore dopo", vtimes(k)
		k=k+1
	end do
end do

!do k=1,9
!	write (*,*), "vettore dopo dopo", vtimes(k)
!enddo

call min_find( vtimes, mtime )

!write (*,*), "Min delle 9 copie: ", mtime
endsubroutine ct_copies

!subroutine evolve( rx, ry, vx, vy, ctime )
!use prec_def
!use constants
!implicit none
!real( long ), dimension( 1:NOD ), intent( in ) :: vx, vy
!real( long ), intent( in ) :: ctime
!real( long ), dimension( 1:NOD ), intent( inout ) :: rx, ry
! local var
!integer k

!do k = 1, NOD
!	rx( k ) = rx( k ) + vx( k ) * ctime
!	ry( k ) = ry( k ) + vy( k ) * ctime
!	! projection inside the box
!	rx( k ) = rx( k ) - floor( rx( k ) )
!	ry( k ) = ry( k ) - floor( ry( k ) )
!end do

!end subroutine evolve

function distance( x1, y1, x2, y2 )
use prec_def
implicit none
real( long ), intent( in ) :: x1, x2, y1, y2
real( long ) :: distance
! local var
real( long ) :: a1, a2, b1, b2
a1 = x1; a2 = x2; b1 = y1; b2 = y2

if( abs( a1-a2 ) > 0.5 ) a2 = a2 - sign( 1.0_long, a2 - a1 )
if( abs( b1-b2 ) > 0.5 ) b2 = b2 - sign( 1.0_long, b2 - b1 )

distance = sqrt( (a2 - a1 )**2 + ( b2 - b1 )**2 )

endfunction distance


!subroutine up_speed( rx, ry, vx, vy, p1, p2 ) ! updates the speeds of the disks after the evolution
!use prec_def
!use constants
!implicit none
!real( long ), dimension( 1:NOD ), intent( in ) :: rx, ry
!real( long ), dimension( 1:NOD ), intent( inout ) :: vx, vy
!integer, intent( in ) :: p1, p2
! local variables
!real(long), dimension( 2 ) :: r1, r2, v1, v2, vers, v12 ! versor r12 = ( r2 - r1 ) / norm
!real(long) :: norm, d, distance ! norm( r2 - r1 )

! index of the vector: i=1 x-component, i=2 y-component
!r1( 1 ) = rx( p1 ); r1( 2 ) = ry( p1 )
!r2( 1 ) = rx( p2 ); r2( 2 ) = ry( p2 )
!v1( 1 ) = vx( p1 ); v1( 2 ) = vy( p1 )
!v2( 1 ) = vx( p2 ); v2( 2 ) = vy( p2 )

! chooses the minimum distance between two points in the box
!if( abs( r2(1) - r1(1) ) > 0.5 ) then
!	r2(1) = r2(1) - sign( 1.0_long, r2(1) - r1(1) )
!endif

!if( abs( r2(2) - r1(2) ) > 0.5 ) then
!	r2(2) = r2(2) - sign( 1.0_long, r2(2) - r1(2) )
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
