program practica5
   implicit none

   interface
      subroutine accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)
         integer, intent(in) :: ndades
         double precision, intent(in) :: xlow, xhigh, cotasup
         double precision, dimension(ndades) :: numeros
         interface
            double precision function funcio(x)
               double precision, intent(in) :: x
            end function funcio
         end interface

      end subroutine accepta

      subroutine histograma(ndat, xdata, xa, xb, nbox, xhis, vhis, errhis, boxsize)
         integer, intent(in) :: ndat
         double precision, intent(in), dimension(ndat) :: xdata
         double precision, intent(in) :: xa, xb
         integer, intent(in) :: nbox
         double precision, intent(out), dimension(nbox) :: xhis, errhis, vhis
         double precision, intent(out) :: boxsize
      end subroutine histograma

      double precision function funcio_p(x)
         double precision, intent(in) :: x
      end function funcio_p

      double precision function funcio_g(x)
         double precision, intent(in) :: x
      end function funcio_g

      subroutine boxmuller(ndades,numeros, mu, sigma)
         integer, intent(in) :: ndades
         double precision, dimension(ndades), intent(out) :: numeros
         double precision, intent(in) :: mu, sigma
      end subroutine boxmuller

      double precision function standard_deviation(values, count)
         integer, intent(in) :: count
         double precision, dimension(count), intent(in) :: values
      end function standard_deviation

      double precision function average(values, count)
         integer, intent(in) :: count
         double precision, dimension(count), intent(in) :: values
      end function average


      double precision function variance(values, count)
         integer, intent(in) :: count
         double precision, dimension(count), intent(in) :: values
      end function variance

   end interface

   double precision, dimension(50000) :: numeros_p
   double precision, dimension(20000) :: numeros_g
   double precision, dimension(100) :: xhis, vhis, errhis
   double precision :: boxsize

   integer :: i
   double precision :: integral, N

   COMMON /CONSTS/PI
   double precision :: PI
   COMMON /FUNC_P/L
   double precision :: L
   COMMON /FUNC_G/SIGMA
   double precision :: SIGMA

   SIGMA = 3
   PI = atan(1d0)*4
   L = 4d0

   call srand(20401021)

   open(1, file="P5-23-24-res.dat")

   ! -- 1 --
   call accepta(50000, numeros_p, -L*PI, L*PI, 1d0/(L*PI), funcio_p)
   call histograma(50000, numeros_p, -L*PI, L*PI, 100, xhis, vhis, errhis, boxsize)

   write(1, *) "# funcio p"
   write(1, *) "# xhis, vhis, errhis"
   do i = 1, 100
      write(1, "(E20.12,E20.12,E20.12)") xhis(i), vhis(i), errhis(i)
   end do
   write(1, *)
   write(1, *)

   call simpsonstresvuit(-L*PI, L*PI/2d0, 12, funcio_p, integral)
   call simpsonstresvuit(-L*PI, L*PI, 12, funcio_p, N)

   write(1, *) "# probabilitat, N"
   write(1, "(E20.12, E20.12)") integral, N
   write(1, *)
   write(1, *)

   ! -- 2 --
   call boxmuller(20000, numeros_g, 0d0, SIGMA)
   call histograma(20000, numeros_g, -4d0*SIGMA, 4d0*SIGMA, 100, xhis, vhis, errhis, boxsize)

   write(1, *) "# funcio g"
   write(1, *) "# xhis, vhis, errhis"
   do i = 1, 100
      write(1, "(E20.12,E20.12,E20.12)") xhis(i), vhis(i), errhis(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# avg, var, desvest"
   write(1, "(E20.12,E20.12,E20.12)") average(numeros_g, 20000), variance(numeros_g, 20000), standard_deviation(numeros_g, 20000)

   close(1)

end program practica5

double precision function funcio_p(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   COMMON /FUNC_P/L
   double precision :: L
   COMMON /CONSTS/PI
   double precision :: PI

   retval = (dsin(x/L)**2)/(L*PI)

end function funcio_p

subroutine boxmuller(ndades,numeros,mu, sigma)
   implicit none

   integer, intent(in) :: ndades
   double precision, dimension(ndades), intent(out) :: numeros
   double precision, intent(in) :: mu, sigma

   COMMON /CONSTS/PI

   double precision :: PI
   double precision :: x1, x2, w1, w2
   integer :: k

   do k = 1, ndades, 2
      w1 = dsqrt(-2d0 * dlog(dble(rand())))
      w2 = 2d0 * PI * dble(rand())

      x1 = mu + sigma * w1 * dcos(w2)
      x2 = mu + sigma * w1 * dsin(w2)

      numeros(k) = x1
      numeros(k+1) = x2
   end do
end subroutine boxmuller

double precision function funcio_g(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   COMMON /FUNC_G/SIGMA
   double precision :: SIGMA

   COMMON /CONSTS/PI
   double precision :: PI

   retval = dexp(-x**2d0/(2d0*SIGMA**2))/dsqrt(2d0*PI*SIGMA**2)
end function funcio_g

subroutine simpsonstresvuit(x1, x2, k, func, integral)
   implicit none

   interface
      double precision function func(x)
         implicit none
         double precision, intent(in) :: x
      end function func
      double precision function simpsonstresvuitfunc(x1, x2, k, func)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function simpsonstresvuitfunc
   end interface

   double precision, intent(in) :: x1, x2
   integer, intent(in) :: k
   double precision, intent(out) :: integral

   integral = simpsonstresvuitfunc(x1, x2, k, func)
end subroutine simpsonstresvuit

double precision function simpsonstresvuitfunc(x1, x2, k, func) result(retval)
   implicit none

   interface
      double precision function func(x)
         implicit none
         double precision, intent(in) :: x
      end function func
      double precision function calc_step(a, b, k1)
         implicit none
         double precision, intent(in) :: a, b
         integer, intent(in) :: k1
      end function calc_step

   end interface

   double precision, intent(in) :: x1, x2
   integer, intent(in) :: k
   integer(kind=8) :: i, intervals
   double precision :: step, integral
   intervals = 3**k
   step = calc_step(x1, x2, k)

   integral = func(x1)

   do i = 1, intervals-3, 3
      integral = integral + 3 * func(x1+step*i) &
         + 3 * func(x1+step*(i+1)) &
         + 2 * func(x1+step*(i+2))
   end do

   integral = integral + 3*(func(x2-2*step) + func(x2-step)) + func(x2)
   retval = 3d0/8d0 * step * integral
end function simpsonstresvuitfunc

double precision function calc_step(a, b, k) result(retval)
   implicit none
   double precision, intent(in) :: a, b
   integer, intent(in) :: k
   retval = (b-a)/dble(3**k)
end function calc_step

subroutine accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)
   implicit none

   integer, intent(in) :: ndades
   double precision, intent(in) :: xlow, xhigh, cotasup
   double precision, dimension(ndades) :: numeros

   interface
      double precision function funcio(x)
         double precision, intent(in) :: x
      end function funcio
      double precision function average(values, count) result(retval)
         integer, intent(in) :: count
         double precision, dimension(count), intent(in) :: values
      end function average
      double precision function variance(values, count) result(retval)
         integer, intent(in) :: count
         double precision, dimension(count), intent(in) :: values
      end function variance
      double precision function standard_deviation(value, normalization_const, weight) result(retval)
         double precision, intent(in) :: value, normalization_const, weight
      end function standard_deviation

   end interface


   double precision :: x1, x2
   integer :: k

   k = 1
   do while (k <= ndades)
      x1 = (xhigh-xlow)*rand() + xlow
      x2 = cotasup * rand()
      if (funcio(x1) > x2) then
         numeros(k) = x1
         k = k + 1
      end if
   end do

end subroutine accepta

subroutine histograma(ndat, xdata, xa, xb, nbox, xhis, vhis, errhis, boxsize)
   implicit none

   interface
      double precision function map(x, a1, b1, a2, b2) result(retval)
         double precision, intent(in) :: x, b1, a1, a2, b2
      end function map
      double precision function standard_deviation(value, total1, weight) result(retval)
         implicit none
         double precision, intent(in) :: value, weight, total1
      end function standard_deviation
   end interface

   integer, intent(in) :: ndat
   double precision, intent(in), dimension(ndat) :: xdata
   double precision, intent(in) :: xa, xb
   integer, intent(in) :: nbox
   double precision, intent(out), dimension(nbox) :: xhis, errhis, vhis
   double precision, intent(out) :: boxsize

   integer :: i, k, total

   boxsize = (xb - xa)/nbox
   k = 1
   vhis = 0
   total = 0
   do i = 1, ndat
      if (xa <= xdata(i) .and. xdata(i) <= xb) then
         k = int(map(xdata(i), xa, xb, 1d0, dble(nbox+1)))
         if (k == nbox +1) then
            k = nbox
         end if
         vhis(k) = vhis(k) + 1
         total = total + 1
      end if
   end do

   do i = 1, nbox
      xhis(i) = xa + (i-0.5d0)*boxsize
      errhis(i) = standard_deviation(vhis(i), dble(total), boxsize) ! subrutina exemple
      vhis(i) = vhis(i)/(boxsize*total) ! normalització
   end do

end subroutine histograma

double precision function deviation(values, count) result(retval)
   implicit none
   interface
      double precision function variance(values1, count1)
         integer, intent(in) :: count1
         double precision, dimension(count1), intent(in) :: values1
      end function variance
   end interface

   integer, intent(in) :: count
   double precision, dimension(count), intent(in) :: values

   retval = dsqrt(variance(values, count))
end function deviation

double precision function variance(values, count) result(retval)
   implicit none

   interface
      double precision function average(values1, count1)
         integer, intent(in) :: count1
         double precision, dimension(count1), intent(in) :: values1
      end function average
   end interface

   integer, intent(in) :: count
   double precision, dimension(count), intent(in) :: values
   double precision, dimension(count) :: distances
   double precision :: avg
   integer :: i

   avg = average(values, count)

   do i = 1, count
      distances(i) = abs(values(i) - avg)**2
   end do

   retval = average(distances, count)
end function variance

double precision function average(values, count) result(retval)
   implicit none
   integer, intent(in) :: count
   double precision, dimension(count), intent(in) :: values
   retval = sum(values, 1) / dble(count)
end function average

double precision function map(x, a1, b1, a2, b2) result(retval)
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map
