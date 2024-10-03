program PRACTICA2
   implicit none

   interface
      function linear_interpolation_series(xs, ys, length, steps_between) result(retval)
         implicit none
         integer, intent(in) :: steps_between, length
         real, dimension(length), intent(in) :: xs, ys
         real, dimension(2, (length-1)*steps_between+1) :: retval

      end function linear_interpolation_series
   end interface

   real, dimension(4) :: x
   real, dimension(500) :: ts, xs_2
   real :: t = 0, step = 0.01, temp

   integer :: i
   real, dimension(2, 1201) :: xy_interpl

   ! Apartat 3
   open(1, file='P2-23-24-res1.dat')
   do while (t <= 5)
      call posipisto(t, x)
      write(1, "(F4.2, F16.6, F16.6, F16.6, F16.6)") t, x(1), x(2), x(3), x(4)
      t = t+step
   end do
   close(1)

   ! Apartat 5
   open(2, file='P2-23-24-res1.dat')
   do i = 1, 500
      read(2, *) ts(i), temp, xs_2(i)
   end do
   close(2)

   xy_interpl = linear_interpolation_series(ts, xs_2, 300, 4)

   open(3, file='P2-23-24-res2.dat')
   do i = 1, (299*4)
      write(3, *) xy_interpl(1, i), xy_interpl(2, i)
   end do
   close(3)
end program PRACTICA2

subroutine posipisto(t, x)
   implicit none
   real, intent(in) :: t
   real, dimension(4), intent(out) :: x
   integer :: i
   real :: calc_posipisto

   do i = 1, 4
      x(i) = calc_posipisto(i, t)
   end do
end subroutine posipisto

real function calc_posipisto(i, t) result(retval)
   integer, intent(in) :: i
   real, intent(in) :: t

   real :: freqmano, radimano, freq, radi
   real :: L = 25., W_0 = 4.8
   freq = freqmano(W_0, i)
   radi = radimano(L, i)

   retval = radi * cos(freq*t)+sqrt(L**2 - radi**2*sin(freq*t)**2)
end function calc_posipisto

real function freqmano(w_0, k) result(retval)
   implicit none
   real, intent(in) :: w_0
   integer, intent(in) :: k

   retval = w_0 * (k/3.25 + 1)
end function freqmano

real function radimano(L, k) result(retval)
   implicit none
   integer, intent(in) :: k
   real, intent(in) :: L

   retval = L - 0.15 - 0.3*(k-1)
end function radimano

!real, dimension(2, length*steps_between)
function linear_interpolation_series(xs, ys, length, steps_between) result(retval)
   implicit none

   interface
      function linear_interpolation(x_0, y_0, x_1, y_1, x_steps) result(r)
         implicit none
         real, intent(in) :: x_0, y_0, x_1, y_1
         integer, intent(in) :: x_steps
         real, dimension(2, x_steps+1) :: r
      end function linear_interpolation
   end interface

   integer, intent(in) :: steps_between, length
   real, dimension(length), intent(in) :: xs, ys
   real, dimension(2, (length-1)*steps_between+1) :: retval

   real, dimension(2, steps_between+1) :: temp
   integer :: i, j


   do i = 1, length-1
      temp = linear_interpolation(xs(i), ys(i), xs(i+1), ys(i+1), steps_between)

      do j = 1, steps_between
         retval(1, (i-1)*steps_between+j) = temp(1, j)
         retval(2, (i-1)*steps_between+j) = temp(2, j)
      end do
   end do

   retval(1, (length-1)*steps_between+1) = xs(length)
   retval(2, (length-1)*steps_between+1) = ys(length)


end function linear_interpolation_series

! real, dimension(x_steps, 2)
function linear_interpolation(x_0, y_0, x_1, y_1, x_steps) result(retval)
   implicit none

   interface
      real function rect(x, m, b) result(y)
         implicit none
         real, intent(in) :: x, m, b
      end function rect
      real function calc_slope(x0, y0, x1, y1) result(m)
         implicit none
         real, intent(in) :: x0, y0, x1, y1
      end function calc_slope
   end interface


   real, intent(in) :: x_0, y_0, x_1, y_1
   integer, intent(in) :: x_steps
   real, dimension(2, x_steps+1) :: retval

   integer :: i
   real :: slope, independent_term ! m, b
   real :: increment, temp_x


   increment = (x_1 - x_0)/x_steps

   slope = calc_slope(x_0, y_0, x_1, y_1)
   independent_term = y_0 - slope*x_0

   retval(1, 1) = x_0
   retval(2, 1) = y_0

   do i = 2, x_steps
      temp_x = x_0 + (i-1)*increment
      retval(1, i) = temp_x
      retval(2, i) = rect(temp_x, slope, independent_term)
   end do

   retval(1, x_steps+1) = x_1
   retval(2, x_steps+1) = y_1

end function linear_interpolation

real function calc_slope(x_0, y_0, x_1, y_1) result(m)
   implicit none
   real, intent(in) :: x_0, y_0, x_1, y_1
   m = (y_1 - y_0)/(x_1 - x_0)
end function calc_slope

real function rect(x, m, b) result(y)
   implicit none
   real, intent(in) :: x, m, b
   y = m*x + b
end function rect
