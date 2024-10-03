program test
   implicit none

   interface
      subroutine euler1(derivative, xini, xend, yini, steps, result_xs, result_ys, result_dys)
         implicit none

         double precision, intent(in) :: xini, xend, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_ys, result_dys

         interface
            double precision function derivative(x, y)
               implicit none
               double precision, intent(in) :: x, y
            end function derivative
         end interface
      end subroutine euler1

      subroutine euler2(double_derivative, xini, xend, yini, iyini, steps, result_xs, result_dys, result_ddys, result_ys)
         implicit none

         double precision, intent(in) :: xini, xend, yini, iyini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ddys, result_ys

         interface
            double precision function double_derivative(x, y)
               implicit none
               double precision, intent(in) :: x, y
            end function double_derivative
         end interface
      end subroutine euler2

      double precision function testxy(x, y) 
          implicit none
          double precision, intent(in) :: x, y
      end function testxy
   end interface

   double precision, dimension(4) :: result_x, result_y, result_dy
   integer :: i

   call euler1(testxy, 0.13d0, 0.14d0, 0.32d0, 4, result_x, result_y, result_dy)
   do i = 1, 4
      write(*, "(E20.12, E20.12, E20.12, E20.12)") result_x(i), result_y(i), result_dy(i)
   end do

end program test

double precision function testxy(x, y) result(retval)
   implicit none
   double precision, intent(in) :: x, y
   retval = dsin(x) - dlog(y)
end function testxy

subroutine euler1(derivative, xini, xend, yini, steps, result_xs, result_ys, result_dys)
   implicit none

   double precision, intent(in) :: xini, xend, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_ys, result_dys

   interface
      double precision function derivative(x, y)
         implicit none
         double precision, intent(in) :: x, y
      end function derivative
   end interface

   integer :: i
   double precision :: h, current_x, current_y

   current_x = xini
   current_y = yini


   h = (xend - xini) / dble(steps)

   do i = 1, steps
      result_dys(i) = derivative(current_x, current_y)
      result_xs(i) = xini + h*i
      result_ys(i) = current_y + h* result_dys(i)
      current_x = result_xs(i)
      current_y = result_ys(i)
   end do

end subroutine euler1

subroutine euler_doble_derivada(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ddys, result_ys)
   implicit none

   double precision, intent(in) :: xini, xend, dyini, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ddys, result_ys

   interface
      double precision function double_derivative(x, y)
         implicit none
         double precision, intent(in) :: x, y
      end function double_derivative
      subroutine euler1(derivative, xini, xend, yini, steps, result_xs, result_ys, result_dys)
         implicit none

         double precision, intent(in) :: xini, xend, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_ys, result_dys

         interface
            double precision function derivative(x, y)
               implicit none
               double precision, intent(in) :: x, y
            end function derivative
         end interface
      end subroutine euler1
   end interface

   integer :: i
   double precision :: h


   h = (xend - xini) / dble(steps)

   call euler1(double_derivative, xini, xend, dyini, steps, result_xs, result_dys, result_ddys)

   result_ys(1) = yini

   do i = 2, steps
      result_ys(i) = result_ys(i-1) + h*result_dys(i)
   end do

end subroutine euler_doble_derivada

