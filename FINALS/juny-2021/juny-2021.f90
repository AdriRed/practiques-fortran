program juny2021
   implicit none

   interface
      double precision function trapezoidalrule(x1, x2, k, func)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function trapezoidalrule
      double precision function euler_mclaurin(x1, x2, k, func, dfunc)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
            double precision function dfunc(x)
               implicit none
               double precision, intent(in) :: x
            end function dfunc
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function euler_mclaurin

      subroutine montecarlosample(func, distr, random_nums1, count, integral, error)
         integer, intent(in) :: count
         double precision, intent(in), dimension(count) :: random_nums1
         double precision, intent(out) :: integral, error
         interface
            double precision function func(x)
               double precision, intent(in) :: x
            end function func
            double precision function distr(x)
               double precision, intent(in) :: x
            end function distr
         end interface
      end subroutine montecarlosample

      double precision function generate_random_distributed_number(distr, a, b, cotasup) result(retval)
         implicit none
         interface
            double precision function distr(x1)
               double precision, intent(in) :: x1
            end function distr
         end interface
         double precision, intent(in) :: a, b, cotasup
      end function generate_random_distributed_number

      double precision function membrane_force(x)
         implicit none
         double precision, intent(in) :: x
      end function membrane_force
      double precision function dmembrane_force(x)
         implicit none
         double precision, intent(in) :: x
      end function dmembrane_force
      double precision function rho(x) result(retval)
         implicit none
           double precision, intent(in) :: x
      end function rho
   end interface

   common /membrane/alpha, k, a
   double precision :: alpha = 1d-3, k = 10d0, a = 40d0
   common /probability/a_prob
   double precision :: a_prob = 1.8d6
   double precision, dimension(100*10**6) :: random_nums

   double precision :: xo = 0, xf = 600, integral1, integral2, master_integral, error
   integer :: m

   open(1, file="resE1.dat")
   write(1, *) "# calcul d'integrals "
   write(1, *) '# m, trapezis, error_tr, mclaurin, error_mc'

   master_integral = euler_mclaurin(xo, xf, 24, membrane_force, dmembrane_force)

   do m = 5, 24
      write(*, *) "Calculating for m = ", m
      integral1 = trapezoidalrule(xo, xf, m, membrane_force)
      integral2 = euler_mclaurin(xo, xf, m, membrane_force, dmembrane_force)
      write(1, "(I2, E20.12, E20.12, E20.12, E20.12)") m, &
         integral1, abs(integral1-master_integral), &
         integral2, abs(integral2-master_integral)
   end do

   close(1)

   open(2, file="resE2.dat")
   write(*, *) "Generating random numbers"
   do m = 1, 100*10**6
      random_nums(m) = generate_random_distributed_number(rho, xo, xf, rho(xf))
   end do
   write(2, *) "# montecarlo"
   write(2, *) "# m, integral"
   ! do m = 1, 100
   !    write(*, *) "Montecarlo for m =", m
   call montecarlosample(membrane_force, rho, random_nums, 100*10**6, integral1, error)
   write(2, "(I3, E20.12)") m, integral1
   ! end do

   close(2)

end program juny2021


double precision function membrane_force(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /membrane/alpha, k, a
   double precision :: alpha, k, a

   retval = k*x/(alpha*dexp(x/a)+1)
end function membrane_force

double precision function dmembrane_force(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /membrane/alpha, k, a
   double precision :: alpha, k, a

   retval = k/(alpha * dexp(x/a)+1d0) - alpha*x*k*dexp(x/a)/(a*(alpha * dexp(x/a)+1d0)**2)
end function dmembrane_force

double precision function rho(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /membrane/alpha, k, a
   double precision :: alpha, k, a
   common /probability/a_prob
   double precision :: a_prob
   retval = k*x/a_prob
end function rho

subroutine montecarlosample(func, distr, random_nums, count, integral, error)
   implicit none
   double precision, intent(in), dimension(count) :: random_nums
   double precision, intent(out) :: integral, error
   integer, intent(in) :: count
   interface
      double precision function func(x1)
         double precision, intent(in) :: x1
      end function func
      double precision function distr(x1)
         double precision, intent(in) :: x1
      end function distr
   end interface

   double precision :: sum, sum2, x, di
   integer :: i

   sum = 0d0
   sum2 = 0d0

   do i = 1, count
      x = random_nums(i)
      sum = sum + func(x)/distr(x)
      sum2 = sum2 + (func(x)**2/distr(x)**2)
      if (mod(i, 10**6) == 0) then
         write(2, "(I3, E20.12)") i/10**6, sum/dble(i)
         write(*, *) "Montecarlo for m =", i/10**6
      end if
   end do

   di = dble(count)

   integral = sum / di
   error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
end subroutine montecarlosample

double precision function trapezoidalrule(x1, x2, k, func) result(retval)
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
   integer :: i, intervals
   double precision :: step, integral
   intervals = 2**k
   step = calc_step(x1, x2, k)

   integral = (func(x2) + func(x1))/2d0

   do i = 1, intervals-1
      integral = integral + func(x1 + step*i)
   end do
   retval = step*integral
end function trapezoidalrule

double precision function euler_mclaurin(x1, x2, k, func, dfunc) result(retval)
   implicit none

   interface
      double precision function func(x)
         implicit none
         double precision, intent(in) :: x
      end function func
      double precision function dfunc(x)
         implicit none
         double precision, intent(in) :: x
      end function dfunc
      double precision function calc_step(a, b, k) result(retval)
         implicit none
         double precision, intent(in) :: a, b
         integer, intent(in) :: k
      end function calc_step
      double precision function trapezoidalrule(x1, x2, k, func) result(retval)
         implicit none

         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface

         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function trapezoidalrule
   end interface

   double precision, intent(in) :: x1, x2
   integer, intent(in) :: k
   double precision :: h

   h = calc_step(x1, x2, k)

   retval = trapezoidalrule(x1, x2, k, func) - h**2/12d0 * (dfunc(x2) - dfunc(x1))

end function euler_mclaurin

double precision function calc_step(a, b, k) result(retval)
   implicit none
   double precision, intent(in) :: a, b
   integer, intent(in) :: k
   retval = (b-a)/dble(2**k)
end function calc_step

double precision function generate_random_distributed_number(distr, a, b, cotasup) result(retval)
   implicit none

   interface
      double precision function distr(x1)
         double precision, intent(in) :: x1
      end function distr
      double precision function map(x1, a1, b1, a2, b2)
         double precision, intent(in) :: x1, a1, b1, a2, b2
      end function map
   end interface

   double precision, intent(in) :: a, b, cotasup



   double precision :: x, y

   x = map(dble(rand()), 0d0, 1d0, a, b)
   y = cotasup * rand()

   do while (y > distr(x))
      x = map(dble(rand()), 0d0, 1d0, a, b)
      y = cotasup * rand()
   end do

   retval = x
end function generate_random_distributed_number

double precision function map(x, a1, b1, a2, b2) result(retval)
   implicit none
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map

