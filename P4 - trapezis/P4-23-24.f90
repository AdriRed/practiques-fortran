program practica2
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
      double precision function high_order_trapezoidalrule(x1, x2, k, func)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function high_order_trapezoidalrule

      double precision function simpsonstresvuit(x1, x2, k, func)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k
      end function simpsonstresvuit
      double precision function ykohoutek(x)
         implicit none
         double precision, intent(in) :: x
      end function ykohoutek
      double precision function calc_step(a, b, k1)
         implicit none
         double precision, intent(in) :: a, b
         integer, intent(in) :: k1
      end function calc_step
   end interface

   common /kohoutek/a, b
   double precision :: a, b
   common /consts/PI
   double precision :: PI

   double precision :: area_exacte, step, x1, x2, area_s, area_t
   integer :: k

   a = 508.633d0 ! x10^6 km
   b = 429.074d0 ! x10^6 km
   PI = atan(1d0)*4

   area_exacte = a*b*(3*sqrt(3d0)+2*PI)/24d0

   x1 = -4d0*a
   x2= -3.5d0*a

   open(1, file="P4-23-24-res.dat")

   write(1, *) "# Calcul d'arees amb diferents metodes i errors"
   write(1, *) "# h, A_T, A_S, Err_T, Err_S"
   do k = 2, 14
      step = calc_step(x1, x2, k)
      area_t = trapezoidalrule(x1, x2, k, ykohoutek)
      area_s = simpsonstresvuit(x1, x2, k, ykohoutek)
      write(1, "(e20.14, 5x,e20.14, 5x,e20.14, 5x,e20.14, 5x,e20.14)") step, &
         area_t, area_s, &
         abs(area_exacte - area_t), abs(area_exacte - area_s)
   end do

   write(1, *)
   write(1, *)

   write(1, *) "# Calcul d'arees amb trapezis d'ordre superior"
   write(1, *) "# h, A, error"
   do k = 2, 13
      step = calc_step(x1, x2, k)
      area_t = high_order_trapezoidalrule(x1, x2, k, ykohoutek)
      write(1, "(e20.14, 5x,e20.14, 5x, e20.14)") step, &
         area_t, abs(area_exacte - area_t)
   end do

   close(1)


end program practica2

double precision function ykohoutek(x) result(retval)
   double precision, intent(in) :: x
   common /kohoutek/a, b
   double precision :: a, b
   retval = b*sqrt(1d0-((x+4d0*a)/a)**2)
end function ykohoutek

double precision function high_order_trapezoidalrule(x1, x2, k, func) result(retval)
   implicit none

   interface
      double precision function func(x)
         implicit none
         double precision, intent(in) :: x
      end function func
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

   end interface

   double precision, intent(in) :: x1, x2
   integer, intent(in) :: k
   retval = (9d0*trapezoidalrule(x1, x2, k+1, func) - trapezoidalrule(x1, x2, k, func))/8d0
end function high_order_trapezoidalrule

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
   intervals = 3**k
   step = calc_step(x1, x2, k)

   integral = (func(x2) + func(x1))/2d0

   do i = 1, intervals-1
      integral = integral + func(x1 + step*i)
   end do
   retval = step*integral
end function trapezoidalrule

double precision function simpsonstresvuit(x1, x2, k, func) result(retval)
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
end function simpsonstresvuit

double precision function calc_step(a, b, k) result(retval)
   implicit none
   double precision, intent(in) :: a, b
   integer, intent(in) :: k
   retval = (b-a)/dble(3**k)
end function calc_step
