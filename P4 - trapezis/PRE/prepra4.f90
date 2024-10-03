program prepra4
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
      double precision function funcioarea(x)
         implicit none
         double precision, intent(in) :: x
      end function funcioarea
      double precision function funciomassa(x)
         implicit none
         double precision, intent(in) :: x
      end function funciomassa
      double precision function funciomassatemps(x)
         implicit none
         double precision, intent(in) :: x
      end function funciomassatemps
      double precision function calc_step(a, b, k1)
         implicit none
         double precision, intent(in) :: a, b
         integer, intent(in) :: k1
      end function calc_step
   end interface

   common /CONSTS/PI
   double precision :: PI
   common /MASSA/L
   double precision :: L = 35.52d-3/2d0

   integer :: k
   double precision :: area_t, area_s, massa_t, massa_s, temp_area, temp_mass
   double precision :: step


   PI = atan(1d0) * 4
   write(*, *) "Area PI", funcioarea(PI)
   ! Apartat 2
   write(*, *) "Calcul d'area"
   area_t = trapezoidalrule(-PI, PI, 13, funcioarea)
   write(*,*) "Trapezoids: ", area_t

   area_s = simpsonstresvuit(-PI, PI, 13, funcioarea)
   write(*,*) "Simpson 3/8: ", area_s

   write(*, *)

   write(*, *) "Calcul de massa"

   massa_t = trapezoidalrule(-L, L, 13, funciomassa)
   write(*,*) "Trapezoids: ", massa_t

   massa_s = simpsonstresvuit(-L, L, 13, funciomassa)
   write(*,*) "Simpson 3/8: ", massa_s

   ! Apartat 3
   open(1, file="area.dat")
   open(2, file="massa.dat")


   do k = 2, 13
      step = calc_step(-PI, PI, k)

      write(1,*) k, step, abs(area_t - trapezoidalrule(-PI, PI, k, funcioarea)), &
         abs(area_s - simpsonstresvuit(-PI, PI, k, funcioarea))

      step = calc_step(-L, L, k)

      write(2,*) k, step, abs(massa_t - trapezoidalrule(-L, L, k, funciomassa)), &
         abs(massa_s - simpsonstresvuit(-L, L, k, funciomassa))
   end do

   close(1)
   close(2)

   ! Apartat 4
   open(3, file="temps.dat")

   do k = 2, 13
      step = calc_step(asin(-PI/L), asin(PI/L), k)
      write(3,*) step, trapezoidalrule(asin(-PI/L), asin(PI/L), k, funciomassatemps), &
         simpsonstresvuit(asin(-PI/L), asin(PI/L), k, funciomassatemps)
   end do


   close(3)
end program prepra4

double precision function funciomassatemps(t) result(retval)
   implicit none
   interface
      double precision function funciomassa(x)
         implicit none
         double precision, intent(in) :: x
      end function funciomassa
      double precision function funciotemps(x)
         implicit none
         double precision, intent(in) :: x
      end function funciotemps
   end interface
   double precision, intent(in) :: t

   retval = funciomassa(funciotemps(t))

end function funciomassatemps

double precision function funcioarea(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   common /CONSTS/PI
   double precision :: PI

   retval = 0.33d0*(cos(x-2d0)*exp(-x**2-sin(x)))**2*sqrt(PI-x)
end function funcioarea

double precision function funciotemps(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /MASSA/L
   double precision :: L
   retval = L * sin(x)
end function funciotemps

double precision function funciomassa(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   double precision :: rho_0

   common /MASSA/L
   double precision :: L

   rho_0 = 0.72d0

   retval = rho_0*sqrt(1-(x/L)**2)*(1-(x/L))*((x/L)**2+(x/L)+1)
end function funciomassa

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
   integer(kind=8) :: i, intervals
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
end function simpsonstresvuit

double precision function calc_step(a, b, k) result(retval)
   implicit none
   double precision, intent(in) :: a, b
   integer, intent(in) :: k
   retval = (b-a)/real(3**k, kind=8)
end function calc_step
