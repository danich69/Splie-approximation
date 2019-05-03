module spline_approximation
	implicit none

    
	contains
	
        subroutine approximation ( x, y, p, result_array, arg_array )       ! main approximation subroutine                   
        
            integer(4) :: i, j, n, k 		                                
            real(8) :: h, t                   
            real(8), dimension ( : ) :: x, y, p, result_array, arg_array                 
            real(8), dimension ( size ( x ) ) :: S, R, BY	                                    
            real(8), dimension ( size ( x ), 0 : 2 ) :: Q, B, B_transp                       
            real(8), dimension ( size ( x ), 0 : 4 ) :: A, Res
            
            n = size ( x )

            call matrix_constructor( x, A, B)                               ! craete needed tridiagonal matrixes
            
            call diagonal_matrix_constructor( Q, p )                        ! create needed diagonal matrix               
            
            B_transp = transpose_3d(B)                                      ! transpose tridiagonal matrix
            BY = 6*matrix_3d_vector_mult(B, y)                              ! multiply tridiagonal matrix by vector of original data
            
            Res = A + 6 * matmul_3d( convertor( matmul_3d(B, Q) ), B_transp )   ! calculate matrix needed for linear system
            
            call fivediagonal_matrix_algorythm(Res, BY, S)                  ! solve linear system using tridiagonal matrix algorithm
            
            
            R = y - matrix_3d_vector_mult ( convertor( matmul_3d( Q, B_transp) ), S)   !create vector needed for calculating spline value
            
                     
            i = 1    
            h = x(i+1) - x(i)
                
                do j = 1, size( arg_array )
                    
                    if( arg_array(j) == x(n) ) then
                        exit
                    endif
                    
                    if( arg_array(j) >= x(i+1) ) then
                    
                        i = i+1
                        h = x(i+1) - x(i)
                        
                    endif
                    
                    t = ( arg_array(j) - x(i) )/h
                    result_array(j) = R(i)*(1-t) + R(i+1)*t - t*(1-t)*(h**2)*((2-t)*S(i)+(1+t)*S(i+1))/6.0      ! calculation of spline value in xi point
                
                end do
                
                result_array(size(result_array)) = R(n)                         ! value of last approximation point
            
		end subroutine
        
        subroutine matrix_constructor(X, A, B)                          ! subroutine which constructs 2 tridiagonal matrixes
    
            real(8) :: x(0:), A(0:,0:), B(0:,0:)
            integer :: n, i
    
            n = size(x) - 1
      
            A=0
            B=0
            
            A(0,2) = 2 * (x(1) - x(0))
            A(0,3) = 0
            A(1,2) = 2 * ( x(2) - x(0) )
           
            A(1,3) = x(2) - x(1)
            B(1,0) = 1/( x(1) - x(0) )
            B(1,1) = -1/( x(1) - x(0) ) + (-1)/( x(2) - x(1) )
            B(1,2) = 1/( x(2) - x(1) )
            
            do i = 2, n-2
        
                A(i,1) = x(i) - x(i-1)
                A(i,2) = 2 * (x(i+1) - x(i-1))
                A(i,3) = x(i+1) - x(i)
                B(i,0) = 1/( x(i) - x(i-1) )
                B(i,1) = -( 1/( x(i) - x(i-1) ) ) - ( 1/( x(i+1) - x(i) ) )
                B(i,2) = 1/( x(i+1) - x(i) )
                
            end do
    
            A(n,2) = 2 * ( x(n) - x(n-1) )
            A(n-1,2) = 2 * ( x(n) - x(n-2) )
            A(n-1,1) = x(n-1) - x(n-2)
            
            B(n-1,0) = 1/( x(n-1) - x(n-2) )
            B(n-1,1) = -(1/( x(n-1) - x(n-2)) ) - ( 1/( x(n) - x(n-1) ) )
            B(n-1,2) = 1/( x(n) - x(n-1) )
            
        end subroutine matrix_constructor
        
        subroutine diagonal_matrix_constructor( Q, p )                          ! subroutine which constructs diagonal matrix
      
            real(8) :: p(1:), Q(:,0:)
            integer :: n, i
            
            n = size(p)
            
            Q = 0
            
            do i = 1, n
            
                Q(i,1) = 1/p(i)
                
            end do
            
        end subroutine diagonal_matrix_constructor
        
        function transpose_3d(A) result(B)                                  ! function which transposes tridiagonal matrix
      
            real(8) :: A(1:,1:), B(1:size(A,1), 1:3)
            integer :: i
        
            B=0
        
            do i = 1, size(A,1) - 1
        
                B(i+1,1) = A(i,3)
                B(i,2) = A(i,2)
                B(i,3) = A(i+1,1)
        
            end do
        
            B(size(A,1),2) = A(size(A,1),2)
        
        end function transpose_3d
        
        function matmul_3d(A, B) result(C)                                  ! function which multipies 2 tridiagonal matrixes
     
            real(8), dimension(1:,1:) :: A, B
            real(8), dimension(1:size(A,1),5) :: C
 
            integer n, i
            
            n = size(A, 1)
            
            C = 0.0
            C(1,3) = A(1,2)*B(1,2)+A(1,3)*B(2,1)
            C(1,4) = A(1,2)*B(1,3)+A(1,3)*B(2,2)
            C(1,5) = A(1,3)*B(2,3)

            C(2,2) = A(2,1)*B(1,2)+A(2,2)*B(2,1)
            C(2,3) = A(2,1)*B(1,3)+A(2,2)*B(2,2)+A(2,3)*B(3,1)
            C(2,4) = A(2,2)*B(2,3)+A(2,3)*B(3,2)
            C(2,5) = A(2,3)*B(3,3)

            do i = 3, n-2
       
                C(i,1) = A(i,1)*B(i-1,1)
                C(i,2) = A(i,1)*B(i-1,2)+A(i,2)*B(i,1)
                C(i,3) = A(i,1)*B(i-1,3)+A(i,2)*B(i,2)+A(i,3)*B(i+1,1)
                C(i,4) = A(i,2)*B(i,3)+A(i,3)*B(i+1,2)
                C(i,5) = A(i,3)*B(i+1,3)
                
            end do

            C(n,1) = A(n,1)*B(n-1,1)
            C(n,2) = A(n,1)*B(n-1,2)+A(n,2)*B(n,1)
            C(n,3) = A(n,1)*B(n-1,3)+A(n,2)*B(n,2)

            C(n-1,1) = A(n-1,1)*B(n-2,1)
            C(n-1,2) = A(n-1,1)*B(n-2,2)+A(n-1,2)*B(n-1,1)
            C(n-1,3) = A(n-1,1)*B(n-2,3)+A(n-1,2)*B(n-1,2)+A(n-1,3)*B(n,1)
            C(n-1,4) = A(n-1,2)*B(n-1,3)+A(n-1,3)*B(n,2)
            
        end function matmul_3d
        
        function matrix_3d_vector_mult(A,B) result(C)                                  ! function which multipies tridiagonal matrix on a vector
        
            real(8) :: A(1:,1:), B(1:), C(size(B))
            integer :: i, n
        
            n = size(B)
            
            C = 0
            C(1) = A(1,2)*B(1) + A(1,3)*B(2)
            C(n) = A(n,1)*B(n-1) + A(n,2)*B(n)
        
            do i = 2, n-1
        
                C(i) = A(i,1)*B(i-1) + A(i,2)*B(i) + A(i,3)*B(i+1)
        
            enddo
        
        end function matrix_3d_vector_mult
        
        function convertor(A) result(B)                                                 ! function which converts fivediagonal matrixe into tridiagonal matrix
      
            real(8) :: A(1:,1:), B(1:size(A,1),1:3)
            integer :: i
     
            do i=1, size(A,1)
      
                B(i,1) = A(i,2)
                B(i,2) = A(i,3)
                B(i,3) = A(i,4)
                
            end do
        
        end function convertor
        
                
        subroutine fivediagonal_matrix_algorythm (Sys, d, x )                          ! subroutine which solves linear system using tridiagonal matrix algorithm
		
			integer :: i
			integer :: msize
            real(8), dimension (:,-2:) :: Sys
            real(8), dimension (:) :: d, x
			real, dimension (-1: size(Sys,1) ) :: beta, alpha, p, q, u, r, a, b, c
            
			msize = size(Sys,1)


            a = 0
            a(1:msize) = Sys(:, 0)
			b = 0
            b(1:msize) = Sys(:, 1)
			c = 0
            c(1:msize) = Sys(:, 2)
			alpha = 0
			beta = 0
			p = 0
			q = 0
			r = 0
            x = 0
			
			do i = 1, msize
            
                beta(i) = b(i-1)-p(i-2)*c(i-2)
                alpha(i) = a(i)-p(i-1)*beta(i)-q(i-2)*c(i-2)
                p(i) = (b(i) - q(i-1)*beta(i))/alpha(i)
                q(i) = c(i)/alpha(i)
                r(i) = (d(i) - r(i-1)*beta(i) - r(i-2)*c(i-2))/alpha(i)
                
			end do
			
            x(msize) = r(msize)
            x(msize - 1) = r(msize-1) - p(msize - 1)*x(msize)
            
			do i = msize-2, 1, -1
            
                x(i) = r(i) - p(i)*x(i+1) - q(i)*x(i+2)
                
			end do
            
		end subroutine fivediagonal_matrix_algorythm
        
end module spline_approximation
