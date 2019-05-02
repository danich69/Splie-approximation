program in_out
	use spline_approximation
	implicit none
    	
!	-----------------------------------------------------
!	|	Program approximates input array function		|
!	|			with cubic splines  					|
!	|													|
!	|	Data stored in file "data.dat"					|
!	|													|
!	| 1st row: # number of intervals of array function	|
!	|		1 coloumn: argument values 					|
!	|		2 coloumn: function values 					|
!	|		3 coloumn: weight of point 					|
!	|													|
!	----------------------------------------------------


	real(8), dimension(:), allocatable :: x, y, p, spline, xi           
    real(8) :: h
	integer(4) :: N, i, j, k
	
	open ( unit=1, file='data.dat' )
	open ( unit=2, file='result.dat' )
	
	read(1,'(2x,I6)') N                                      			! read number of intervals
	
	allocate( x(N+1), y(N+1), p(N+1), spline(N*100+1), xi(N*100+1) )
	
	do i=1,N+1
    
		read(1,*) x(i), y(i), p(i)                                      ! read data
        
    end do
    
    h = ( x(N+1) - x(1) )/(100*N)
       
    do i = 0, 100*N                                                     ! create approximation net
        
        xi(i+1) = x(1) + h*i
		
	end do 
    
	call approximation( x, y, p, spline, xi )                           ! call approximation subroutine
    
    
    do i = 1, 100*N+1
        
       write(2,*) xi(i), spline(i)                                      ! print results to result.dat
		
	end do
	
    
	deallocate(x, y, p)
    deallocate(xi, spline)
	
	close(1)

end program in_out
