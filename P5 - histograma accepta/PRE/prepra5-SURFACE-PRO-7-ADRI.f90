program name
   implicit none

   interface
      subroutine histograma(ndat, xdata, xa, xb, nbox, xhis, vhis, errhis, boxsize)
         integer, intent(in) :: ndat
         double precision, intent(in), dimension(ndat) :: xdata
         double precision, intent(in) :: xa, xb
         integer, intent(in) :: nbox
         double precision, intent(out), dimension(nbox) :: xhis, errhis, vhis
         double precision, intent(out) :: boxsize
      end subroutine histograma
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
      double precision function p(x) result(retval)
         double precision, intent(in) :: x
      end function p
   end interface

   double precision, dimension(40000) :: numeros
   double precision, dimension(50) :: xhis, vhis, errhis
   double precision :: boxsize

   integer :: i 

   COMMON /CONSTS/PI
   double precision :: PI

   call srand(20401021)
   open(1, file="P5-23-24.dat")

   PI = atan(1d0)*4

   call accepta(40000, numeros, -PI, PI, 1d0, p)

   call histograma(40000, numeros, -PI, PI, 50, xhis, vhis, errhis, boxsize)

   write(1, *) "# histograma"
   write(1, *) "# xhis, vhis, err"

   do i = 1, 50
      write(1, "(E20.12, E20.12, E20.12)") xhis(i), vhis(i), errhis(i)
   end do

   close(1)
end program name

subroutine sexponencial()
end subroutine sexponencial

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

   ! este codigo no deberia de ir aqui
   write(1, *) "# estadistiques"
   write(1, *) "# avg, var, std.dev"
   write(1, "(E20.12, E20.12, E20.12)") average(numeros, ndades), variance(numeros, ndades), dsqrt(variance(numeros, ndades))
   write(1, *) 
   write(1, *) 


end subroutine accepta

double precision function p(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   COMMON /CONSTS/PI
   double precision :: PI
   retval = 125d0*dexp(PI)*(x*dsin(x))**2*dexp(-dabs(x))/(4d0*(68d0*dexp(PI)-70d0*PI-25d0*PI**2-68d0))
end function p

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

double precision function standard_deviation(value, normalization_const, weight) result(retval)
   implicit none
   double precision, intent(in) :: value, normalization_const, weight
   retval = dsqrt(value/normalization_const * (1d0-value/ normalization_const) ) / (weight * dsqrt(normalization_const))
end function standard_deviation

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

function apply(values, count, fun) result(retval)
   implicit none
   integer, intent(in) :: count
   double precision, dimension(count), intent(in) :: values
   double precision, dimension(count):: retval

   interface
      double precision function fun(x)
         double precision, intent(in) :: x
      end function fun
   end interface

   integer :: i

   do i = 1, count
      retval(i) = fun(values(i))
   end do
end function apply


double precision function map(x, a1, b1, a2, b2) result(retval)
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map


