program practica4
   implicit none

   interface
      subroutine comet_radius(x, f_x, df_x)
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine comet_radius
      subroutine comet_x(x, f_x, df_x)
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine comet_x
      subroutine comet_y(x, f_x, df_x)
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine comet_y
      subroutine function2(x2, f_x2, df_x2)
         double precision, intent(in) :: x2
         double precision, intent(out) :: f_x2, df_x2
      end subroutine function2
      subroutine anomal_eccentricity(x3, f_x, df_x)
         double precision, intent(in) :: x3
         double precision, intent(out) :: f_x, df_x
      end subroutine anomal_eccentricity
      subroutine biseccio(fun, a, b, preci, nitera, valorarrel)
         implicit none
         interface
            subroutine fun(x, f_x, df_x)
               implicit none
               double precision, intent(in) :: x
               double precision, intent(out) :: f_x, df_x
            end subroutine fun
         end interface

         double precision, intent(in) :: a, b, preci
         integer, intent(out) :: nitera
         double precision, intent(out) :: valorarrel
      end subroutine biseccio
      subroutine derivataula(length, xvalues, f_xvalues, df_xvalues)
         integer, intent(in) :: length
         double precision, dimension(length), intent(in) :: xvalues, f_xvalues
         double precision, dimension(length), intent(out) :: df_xvalues
      end subroutine derivataula
      double precision function calc_triangle_area_from_points(x1, y1, x2, y2, x3, y3) result(retval)
         double precision, intent(in) :: x1, y1, x2, y2, x3, y3
      end function calc_triangle_area_from_points
   end interface

   COMMON /COMET/MAJOR_SEMIAXIS, ECCENTRICITY
   double precision :: MAJOR_SEMIAXIS = 17.857619, ECCENTRICITY = 0.967990
   COMMON /AP3/T_H, t
   double precision :: T_H = 75.3
   COMMON /CONSTS/PI
   double precision :: PI = 3.14159265359

   double precision :: step, E, fE, fE2, dfE, dfE2, arrel, t, E_0, dx, dy, area, xfocus, yfocus
   double precision, dimension(80) :: xs, ys
   double precision, dimension(100) :: Es, radii_E, dradii_E
   integer :: i, iteracions

   step = 0.02*PI

   do i = 1, 100
      E = step * (i-1)
      call comet_radius(E, fE, dfE)
      Es(i) = E
      radii_E(i) = fE
   end do

   call derivataula(100, Es, radii_E, dradii_E)


   open(1, file="P3-23-24-res.dat")

   ! Apartat 1
   write(1, *) "# Càlcul del radi"
   do i = 1, 100
      write(1, "(E20.12, E20.12, E20.12)") Es(i), radii_E(i), dradii_E(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat 2
   write(1, *) "# Càlcul d'arrels amb bisecció"

   E = 0.2
   step = 0.1
   do while (E <= 5.8d0)
      call function2(E-step, fE, dfE)
      call function2(E, fE2, dfE2)

      if (fE*fE2 < 0) then
         call biseccio(function2, E-step, E, 1.d-12, iteracions, arrel)
         call comet_radius(arrel, fE, dfE)
         write(1, "(E20.12, E20.12)") arrel, fE
      end if

      E = E+step
   end do

   write(1, *)
   write(1, *)

   ! Apartat 3
   write(1, *) "# Càlcul d'anomalies i posicions"

   t = 0
   E_0 = PI/4.
   step = T_H / 80.
   i = 1
   do while (t <= T_H)
      call newtonraphson(anomal_eccentricity, E_0, 1.d-12, iteracions, E)

      call comet_x(E, xs(i), dx)
      call comet_y(E, ys(i), dy)

      write(1, "(E20.12, E20.12, E20.12, E20.12)") t, E, xs(i), ys(i)

      t = t + step
      i = i + 1
   end do

   xfocus = -(2*MAJOR_SEMIAXIS-sqrt(xs(i)**2+ys(i)**2))
   yfocus = 0d0 ! Està a y = 0
   write(1, "(E20.12, E20.12, E20.12, E20.12)") t, E, xfocus, yfocus
   write(1, *)
   write(1, *)

   ! Apartat extra
   write(1, *) "# Càlcul d'àrea"


   t = 0
   do i = 1, 79
      area = calc_triangle_area_from_points(xs(i), ys(i), xs(i+1), ys(i+1), xfocus, yfocus)
      write(1,  "(E20.12, E20.12)") t, area
      t = t+step
   end do

   write(1, *)
   write(1, *)

end program practica4

subroutine function2(x, f_x, df_x)
   implicit none
   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   COMMON /COMET/MAJOR_SEMIAXIS, ECCENTRICITY
   double precision MAJOR_SEMIAXIS, ECCENTRICITY

   f_x = sin(2*x) * (1-ECCENTRICITY**2) - (cos(x)*(2-ECCENTRICITY**2)-ECCENTRICITY)*sin(x)
   df_x = 2*(1-ECCENTRICITY**2)*cos(2*x) - (cos(x)*(2-ECCENTRICITY**2)-ECCENTRICITY)*cos(x) + sin(x)*(2-ECCENTRICITY**2)*cos(x)
end subroutine function2

subroutine comet_x(x, f_x, df_x)
   implicit none
   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   COMMON /COMET/MAJOR_SEMIAXIS, ECCENTRICITY
   double precision MAJOR_SEMIAXIS, ECCENTRICITY

   f_x = MAJOR_SEMIAXIS*(cos(x)-ECCENTRICITY)
   df_x = -MAJOR_SEMIAXIS*sin(x)
end subroutine comet_x

subroutine anomal_eccentricity(x, f_x, df_x)
   implicit none
   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   COMMON /COMET/MAJOR_SEMIAXIS, ECCENTRICITY
   double precision MAJOR_SEMIAXIS, ECCENTRICITY
   COMMON /AP3/T_H, t
   double precision :: T_H, t
   COMMON /CONSTS/PI
   double precision :: PI

   f_x = 2*PI/T_H*t - ECCENTRICITY*sin(x)-x
   df_x = -ECCENTRICITY*cos(x)-1
end subroutine anomal_eccentricity

subroutine comet_y(x, f_x, df_x)
   implicit none
   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   COMMON /COMET/MAJOR_SEMIAXIS, ECCENTRICITY
   double precision MAJOR_SEMIAXIS, ECCENTRICITY

   f_x = MAJOR_SEMIAXIS*sqrt(1-ECCENTRICITY**2)*sin(x)
   df_x = MAJOR_SEMIAXIS*sqrt(1-ECCENTRICITY**2)*cos(x)
end subroutine comet_y

subroutine comet_radius(x, f_x, df_x)
   implicit none

   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x
   double precision :: xcomet, dxcomet, ycomet, dycomet

   interface
      subroutine comet_x(x1, f_x1, df_x1)
         double precision, intent(in) :: x1
         double precision, intent(out) :: f_x1, df_x1
      end subroutine comet_x
      subroutine comet_y(x2, f_x2, df_x2)
         double precision, intent(in) :: x2
         double precision, intent(out) :: f_x2, df_x2
      end subroutine comet_y
   end interface

   call comet_x(x, xcomet, dxcomet)
   call comet_y(x, ycomet, dycomet)

   f_x = sqrt(xcomet**2+ycomet**2)
   df_x = (1d0/sqrt(xcomet**2+ycomet**2))*(2*xcomet*dxcomet+2*ycomet*dycomet)

end subroutine comet_radius

subroutine mult_funcs(x, f, g, fg_x, dfg_x)
   implicit none
   interface
      subroutine f(x1, f_x1, df_x1)
         double precision, intent(in) :: x1
         double precision, intent(out) :: f_x1, df_x1
      end subroutine f
      subroutine g(x2, f_x2, df_x2)
         double precision, intent(in) :: x2
         double precision, intent(out) :: f_x2, df_x2
      end subroutine g
   end interface

   double precision, intent(in) :: x
   double precision, intent(out) :: fg_x, dfg_x
   double precision :: f_x, df_x, g_x, dg_x

   call f(x, f_x, df_x)
   call g(x, g_x, dg_x)

   fg_x = f_x*g_x
   dfg_x = df_x * g_x + f_x * dg_x

end subroutine mult_funcs

subroutine biseccio(fun, a, b, preci, nitera, valorarrel)
   implicit none

   interface
      subroutine fun(x, f_x, df_x)
         implicit none
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine fun
   end interface

   double precision, intent(in) :: a, b, preci
   integer, intent(out) :: nitera
   double precision, intent(out) :: valorarrel
   double precision :: inner_a, inner_b, f_a, df_a
   double precision :: c, f_c, df_c
   integer :: i

   inner_a = a
   inner_b = b
   i = 0
   do while (abs(inner_b-inner_a) > preci)
      c = (inner_a+inner_b)/2d0
      call fun(c, f_c, df_c)
      if (f_c == 0.) then
         error stop
      end if

      call fun(inner_a, f_a, df_a)
      if (f_c * f_a < 0d0) then
         inner_b = c
      else
         inner_a = c
      endif
      i = i + 1
   end do
   nitera = i
   valorarrel = (inner_a+inner_b)/2d0
end subroutine biseccio

subroutine newtonraphson(fun, xini, preci, nitera, valorarrel)
   implicit none

   interface
      subroutine fun(x, f_x, df_x)
         implicit none
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine fun
   end interface

   double precision, intent(in) :: xini, preci
   integer, intent(out) :: nitera
   double precision, intent(out) :: valorarrel
   double precision :: delta, f_x0, df_x0, x1

   valorarrel = xini
   nitera = 0
   delta = preci +1

   do while (delta > preci)
      call fun(valorarrel, f_x0, df_x0)
      x1 = valorarrel - f_x0/df_x0
      delta = abs(f_x0/df_x0)
      valorarrel = x1
      nitera = nitera + 1
   end do

   valorarrel = valorarrel
end subroutine newtonraphson

subroutine derivataula(length, xvalues, f_xvalues, df_xvalues)
   integer, intent(in) :: length
   double precision, dimension(length), intent(in) :: xvalues, f_xvalues
   double precision, dimension(length), intent(out) :: df_xvalues
   integer :: i
   double precision :: h

   h = xvalues(2) - xvalues(1)

   df_xvalues(1) = (f_xvalues(2) - f_xvalues(1))/h

   do i = 2, length-1
      df_xvalues(i) = (f_xvalues(i+1) - f_xvalues(i-1))/(2*h)
   end do

   df_xvalues(length) = (f_xvalues(length) - f_xvalues(length-1))/h
end subroutine derivataula

double precision function calc_triangle_area_from_points(x1, y1, x2, y2, x3, y3) result(retval)
   double precision, intent(in) :: x1, y1, x2, y2, x3, y3
   double precision :: xmiddle, ymiddle, h, b

   xmiddle = (x2-x1)/2d0
   ymiddle = (y2-y1)/2d0

   b = sqrt((x2-x1)**2+(y2-y1)**2)


   h = sqrt((xmiddle-x3)**2+(ymiddle-y3)**2)

   retval = b*h/2d0

end function calc_triangle_area_from_points

double precision function rect(x, m, b) result(y)
   implicit none
   double precision, intent(in) :: x, m, b
   y = m*x + b
end function rect
