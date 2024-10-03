program PRACTICA2
   implicit none

   COMMON /CONSTS/PI
   COMMON /POSIS/XI,TI

   real, dimension(5) :: x
   real :: t = 0, step = 0.1, temp, interpolated_x, xinterpo, xinterpo0
   real :: PI = 3.14159265359
   real, dimension(81) :: XI, TI

   integer :: i

   ! Apartat 3
   open(1, file='P2-23-24-res1-c.dat')
   do while (t <= 8.)
      call posipisto(t, x)
      write(1, "(F4.2, F16.6, F16.6, F16.6, F16.6, F16.6)") t, x(1), x(2), x(3), x(4), x(5)
      t = t+step
   end do
   close(1)

   ! Apartat 5
   open(2, file='P2-23-24-res1-c.dat')
   do i = 1, 81
      read(2, *) TI(i), temp, temp, temp, XI(i)
   end do
   close(2)
   
   ! Aparatat 6
   step = 0.003
   t = 0.001
   open(3, file='P2-23-24-res2-c.dat')
   do while (t < 6.)
      interpolated_x = xinterpo(t)
      write(3, "(F5.3, F16.6)") t, interpolated_x
      t = t + step
   end do
   close(3)

   t = 0.001
   open(4, file='P2-23-24-res3-c.dat')
   do while (t < 6.)
      interpolated_x = xinterpo0(t)
      write(4, "(F5.3, F16.6)") t, interpolated_x
      t = t + step
   end do
   close(4)

end program PRACTICA2

subroutine posipisto(t, x)
   implicit none
   real, intent(in) :: t
   real, dimension(5), intent(out) :: x
   integer :: i
   real :: calc_posipisto

   do i = 1, 5
      x(i) = calc_posipisto(i, t)
   end do
end subroutine posipisto

real function calc_posipisto(i, t) result(retval)
   integer, intent(in) :: i
   real, intent(in) :: t

   real :: radimano, phi, radi, fase
   real :: L = 18.5, freq = 5.
   radi = radimano(L, i)
   fase = phi(i)

   retval = radi * cos(freq*t+fase)+sqrt(L**2 - radi**2*sin(freq*t + fase)**2)
end function calc_posipisto

real function phi(i) result(retval)
   COMMON /CONSTS/PI

   integer, intent(in) :: i
   real :: PI

   retval = (i/5.)**2*PI
end function phi

real function xinterpo(t) result(retval)
   implicit none
   COMMON /POSIS/XI, TI

   real, intent(in) :: t
   real, dimension(81) :: XI, TI
   integer :: lower_index, find_index_of_lower_x
   real :: slope, independent_term, calc_slope, rect
   real :: x0, y0, x1, y1

   lower_index = find_index_of_lower_x(TI, t)

   x0 = TI(lower_index)
   x1 = TI(lower_index+1)
   y0 = XI(lower_index)
   y1 = XI(lower_index+1)


   slope = calc_slope(x0, y0, x1, y1)
   independent_term = y0 - slope*x0

   retval = rect(t, slope, independent_term)

end function xinterpo

real function xinterpo0(t) result(retval)
   implicit none
   COMMON /POSIS/XI, TI

   real, intent(in) :: t
   real, dimension(81) :: XI, TI
   integer :: lower_index, find_index_of_lower_x
   real :: y0, y1


   lower_index = find_index_of_lower_x(TI, t)

   y0 = XI(lower_index)
   y1 = XI(lower_index+1)


   retval = (y1+y0)/2

end function xinterpo0

real function radimano(L, k) result(retval)
   implicit none
   integer, intent(in) :: k
   real, intent(in) :: L

   retval = L/k - 0.5

end function radimano

integer function find_index_of_lower_x(x_set, x_0) result(retval)
   implicit none
   real, intent(in) :: x_0
   real, dimension(81), intent(in) :: x_set
   real :: current_x
   integer :: i

   do i = 1, 81
      current_x = x_set(i)
      if (current_x > x_0) then
         retval = i-1
         exit
      end if
   end do

end function find_index_of_lower_x

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

 ! -------------- FUNCIONS QUE HE UTILITZAT A LA PREPRA I A LA PRACTICA NO --------------

 ! NOMES QUAN TENIM UNA QUANTIAT ENTERA DE PUNTS ENTRE CADA SAMPLE
 ! NO L'UTILITZO A LA PRACTICA JA QUE TENIM 33.3 punts entre cada sample
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

 ! NOMES QUAN TENIM UNA QUANTIAT ENTERA DE PUNTS ENTRE CADA SAMPLE
 ! NO L'UTILITZO A LA PRACTICA
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
