program prepra9
   implicit none

   interface
      subroutine poisson_equation_algorithm(method_used, width, height, input_matrix, epsilon, rho, h, write_steps, output_matrix)
         implicit none

         interface
            double precision function rho(x, y)
               integer, intent(in) :: x, y
            end function
            subroutine method_used(width, height, rho, h, input_matrix, output_matrix)
               implicit none
               interface
                  double precision function rho(x, y)
                     integer, intent(in) :: x, y
                  end function
               end interface
               integer, intent(in) :: width, height
               double precision :: h
               double precision, intent(in), dimension(width, height) :: input_matrix
               double precision, intent(out), dimension(width, height) :: output_matrix
            end subroutine method_used
         end interface

         integer, intent(in) :: width, height
         double precision, intent(in), dimension(width, height) :: input_matrix
         double precision, intent(in) :: epsilon, h
         double precision, intent(out), dimension(width, height) :: output_matrix
         logical, intent(in) :: write_steps
      end subroutine poisson_equation_algorithm
      double precision function stoves(x, y) result(retval)
         implicit none
         integer, intent(in) :: x, y
      end function stoves
      subroutine jacobi_method(width, height, rho, h, input_matrix, output_matrix)
         implicit none

         interface
            double precision function rho(x, y)
               integer, intent(in) :: x, y
            end function
         end interface

         integer, intent(in) :: width, height
         double precision :: h
         double precision, intent(in), dimension(width, height) :: input_matrix
         double precision, intent(out), dimension(width, height) :: output_matrix
      end subroutine jacobi_method
      subroutine successive_overrelaxation_method(width, height, rho, h, input_matrix, output_matrix)
         implicit none

         interface
            double precision function rho(x, y)
               integer, intent(in) :: x, y
            end function
         end interface

         integer, intent(in) :: width, height
         double precision :: h
         double precision, intent(in), dimension(width, height) :: input_matrix
         double precision, intent(out), dimension(width, height) :: output_matrix
      end subroutine successive_overrelaxation_method
      function initialize_heat_input_matrix(width, height, interior, top, bottom, left, right) result(retval)
         integer, intent(in) :: width, height
         double precision, intent(in) :: interior, top, bottom, left, right
         double precision, dimension(width, height) :: retval
      end function
   end interface

   common /stoves_consts/circular_heat, rectangular_heat, h
   common /overrelaxation/correction
   double precision :: circular_heat, rectangular_heat, correction, Lx, Ly, h, eps
   double precision, dimension(38, 63) :: input, output
   double precision, dimension(3) :: interior_temps
   integer :: max_index_x, max_index_y, i, j
   h = 0.5d0
   circular_heat = 0.8d0
   rectangular_heat = 1.1d0
   correction = 1.52d0
   Lx = 18.5d0
   Ly = 31d0
   max_index_x = 38
   max_index_y = 63
   eps = 1d-4

   interior_temps = (/2d0, 14d0, 230d0/)
   open(1, file="data.dat")
   open(2, file="animation.dat")
   write(2, *) "# animation gif"
   write(2, *) "# x, y, t"
   do i = 1, 3
      input = initialize_heat_input_matrix(max_index_x, max_index_y, interior_temps(i), 4.16d0, 13.15d0, 2.2d0, 35d0)

      write(1, "(A, F5.1)") "# temperature point for T = ", interior_temps(i)
      write(1, *) "# jacobi"
      write(*, *) "Poisson with jacobi START"
      call poisson_equation_algorithm(jacobi_method, max_index_x, max_index_y, input, eps, stoves, h, .false., output)
      write(*, *) "Poisson with jacobi OK"
      write(1, *)
      write(1, *)

      write(1, *) "# successive overrelaxation"
      write(*, *) "Poisson with successive overrelaxation START"
      call poisson_equation_algorithm(successive_overrelaxation_method, max_index_x, max_index_y, input, eps, stoves, h, i == 1, &
         output)
      write(*, *) "Poisson with successive overrelaxation OK"
      write(1, *)
      write(1, *)

   end do

   write(1, *) "# temperature matrix"
   write(1, *) "# x, y, T"

   do i = 1, max_index_x
      do j = 1, max_index_y
         write(1, "(E20.12, E20.12, E20.12)") (i-1)*h, (j-1)*h, output(i, j)
      end do
   end do

   close(1)
   close(2)
end program prepra9

function initialize_heat_input_matrix(width, height, interior, top, bottom, left, right) result(retval)
   implicit none

   integer, intent(in) :: width, height
   double precision, intent(in) :: interior, top, bottom, left, right
   double precision, dimension(width, height) :: retval
   retval(2:width-1, 2:height-1) = interior
   ! T(0, y)
   retval(1, 1:height) = left
   ! T(Lx, y)
   retval(width, 1:height) = right
   ! T(x, 0)
   retval(1:width, 1) = bottom
   ! T(x, Ly)
   retval(1:width, height) = top

end function initialize_heat_input_matrix

subroutine poisson_equation_algorithm(method_used, width, height, input_matrix, epsilon, rho, h, write_steps, output_matrix)
   implicit none

   interface
      double precision function rho(x, y)
         integer, intent(in) :: x, y
      end function
      subroutine method_used(width, height, rho, h, input_matrix, output_matrix)
         implicit none
         interface
            double precision function rho(x, y)
               integer, intent(in) :: x, y
            end function
         end interface
         integer, intent(in) :: width, height
         double precision :: h
         double precision, intent(in), dimension(width, height) :: input_matrix
         double precision, intent(out), dimension(width, height) :: output_matrix
      end subroutine method_used
   end interface

   integer, intent(in) :: width, height
   double precision, intent(in), dimension(width, height) :: input_matrix
   double precision, intent(in) :: epsilon, h
   logical, intent(in) :: write_steps
   double precision, intent(out), dimension(width, height) :: output_matrix

   double precision, dimension(width, height) :: current_matrix
   integer :: i, x, y

   call method_used(width, height, rho, h, input_matrix, output_matrix)
   output_matrix = 0
   i = 1

   current_matrix = input_matrix
   call method_used(width, height, rho, h, current_matrix, output_matrix)

   write(1, "(I5, E20.12)") i, output_matrix(15, 27)
   if (write_steps) then
      do x = 1, width
         do y = 1, height
            write(2, "(E20.12, E20.12, E20.12)") (x-1)*h, (y-1)*h, output_matrix(x, y)
         end do
      end do
      write(2, *)
      write(2, *)
   end if

   do while (maxval(dabs(current_matrix - output_matrix)) > epsilon)
      i = i+1
      current_matrix = output_matrix
      call method_used(width, height, rho, h, current_matrix, output_matrix)
      write(1, "(I5, E20.12)") i, output_matrix(15, 27)
      if (write_steps) then
         do x = 1, width
            do y = 1, height
               write(2, "(E20.12, E20.12, E20.12)") (x-1)*h, (y-1)*h, output_matrix(x, y)
            end do
         end do
         write(2, *)
         write(2, *)   
      end if
   end do
end subroutine poisson_equation_algorithm

subroutine jacobi_method(width, height, rho, h, input_matrix, output_matrix)
   implicit none

   interface
      double precision function rho(x, y)
         integer, intent(in) :: x, y
      end function
   end interface

   integer, intent(in) :: width, height
   double precision :: h
   double precision, intent(in), dimension(width, height) :: input_matrix
   double precision, intent(out), dimension(width, height) :: output_matrix

   integer :: i, j
   output_matrix = input_matrix

   do i = 2, width-1
      do j = 2, height-1
         output_matrix(i, j) = 0.25d0 * (input_matrix(i+1, j) + input_matrix(i-1, j) + &
            input_matrix(i, j+1) + input_matrix(i, j-1) + &
            h**2*rho(i, j))
      end do
   end do
end subroutine jacobi_method

subroutine successive_overrelaxation_method(width, height, rho, h, input_matrix, output_matrix)
   implicit none
   interface
      double precision function rho(x, y)
         integer, intent(in) :: x, y
      end function
   end interface
   integer, intent(in) :: width, height
   double precision :: h
   double precision, intent(in), dimension(width, height) :: input_matrix
   double precision, intent(out), dimension(width, height) :: output_matrix

   common /overrelaxation/correction
   double precision :: correction

   integer :: i, j
   output_matrix = input_matrix
   do i = 2, width-1
      do j = 2, height-1
         output_matrix(i, j) = input_matrix(i, j) + correction * (0.25d0 * (input_matrix(i+1, j) + &
            input_matrix(i, j+1) + output_matrix(i, j-1) + output_matrix(i-1, j) + &
            h**2*rho(i, j)) - input_matrix(i, j))
      end do
   end do

end subroutine successive_overrelaxation_method

double precision function stoves(ix, jy) result(retval)
   implicit none

   integer, intent(in) :: ix, jy

   common /stoves_consts/circular_heat, rectangular_heat, h
   double precision :: circular_heat, rectangular_heat, x, y, h, r

   x = (ix-1) * h
   y = (jy-1) * h
   r = dsqrt((x-9.5d0)**2 + (y-12d0)**2)

   retval = circular_heat*dexp(-(r-2d0)**2/0.35**2)

   if ((12d0 <= x .and. x <= 16d0) .and. (20d0 <= y .and. y <= 28d0)) then
      retval = retval + rectangular_heat
   end if

end function stoves
