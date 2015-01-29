! modules
module prec_def
implicit none
integer, parameter :: long = selected_real_kind( 18,307 )
end module prec_def

program bin
use prec_def

implicit none
real( long ) :: mean, var, mean_m, temp, unc
real( long ), dimension( : ), allocatable :: input
integer :: i,j,k,m, istat, leng, n_bin
integer :: maxb = 100 ! max dimension of bin

print *, "started"
open( unit = 10, file = "pressure.dat", status = "old", access = "sequential", position = "rewind", action = "read" )
open( unit = 100, file = "bin.dat", access = "sequential", position = "rewind", action = "write" )

k = 0
do 
	read( unit = 10, fmt = *, iostat = istat )
	if( istat < 0 ) exit
	k = k + 1
end do
leng = k ! length of the vector containing pressures
print *, leng
allocate( input( 1:leng ) )

close( unit = 10 )
open( unit = 10, file = "pressure.dat", status = "old", access = "sequential", position = "rewind", action = "read" )

k = 1
do 
	read( unit = 10, fmt = 100, iostat = istat ) temp
		100 format( f20.15 )
	if( istat < 0 ) exit
	input( k ) = temp
	k = k + 1
end do


do m = 1, maxb
	n_bin = leng/m  ! number of bins
	var = 0
	mean = 0 ! mean of pressures
	
	
	! calcolo la media sui dati
	do k=1, n_bin * m
		mean = ( mean + input( k ) ) / ( n_bin * m )
	enddo
	
	do i = 1, n_bin * m, m
	mean_m = 0
		do j = i, i -1 + m
			mean_m = ( mean_m + input( j ) ) / m
			! print *, "running"
		end do
		var = ( mean_m - mean )**2 / ( n_bin -1 )
	end do
	
	unc = sqrt( var / ( n_bin ) )
	write( unit = 100, fmt = 101 ) m, unc
		101 format( i20, "	", f20.15 )
end do
print *, "ok"

deallocate( input )

close( unit = 10 )
close( unit = 100 )


print *, "ok"

end program bin
