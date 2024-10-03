program name
   implicit none

   interface
      subroutine euler1(derivative, xini, xend, yini, steps, result_xs, result_ys)
         implicit none

         double precision, intent(in) :: xini, xend, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_ys

         interface
            double precision function derivative(x, y)
               implicit none
               double precision, intent(in) :: x, y
            end function derivative
         end interface
      end subroutine euler1
      subroutine euler2(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
         implicit none

         double precision, intent(in) :: xini, xend, dyini, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ys

         interface
            double precision function double_derivative(x, y, dy)
               implicit none
               double precision, intent(in) :: x, y, dy
            end function double_derivative
         end interface
      end subroutine euler2
      subroutine adams_bashforth1(derivative, xini, xend, yini, steps, result_xs, result_ys)
         implicit none

         double precision, intent(in) :: xini, xend, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_ys

         interface
            double precision function derivative(x, y)
               implicit none
               double precision, intent(in) :: x, y
            end function derivative
         end interface
      end subroutine adams_bashforth1
      subroutine adams_bashforth2(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
         implicit none

         double precision, intent(in) :: xini, xend, dyini, yini
         integer, intent(in) :: steps
         double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ys

         interface
            double precision function double_derivative(x, y, dy)
               implicit none
               double precision, intent(in) :: x, y, dy
            end function double_derivative
         end interface
      end subroutine adams_bashforth2


      double precision function pendol(x, y, dy)
         implicit none
         double precision, intent(in) :: x, y, dy
      end function pendol
      double precision function pendol_simplificat(x, y, dy)
         implicit none
         double precision, intent(in) :: x, y, dy
      end function pendol_simplificat
      double precision function Epoten(x)
         implicit none
         double precision, intent(in) :: x
      end function Epoten
      double precision function Ecine(x)
         implicit none
         double precision, intent(in) :: x
      end function Ecine
   end interface

   double precision :: T_N, w_N
   double precision, dimension(1300) :: pendol_s_t, pendol_s_angle, pendol_s_velocitat
   double precision, dimension(1800) :: pendol_t, pendol_angle, pendol_velocitat
   double precision, dimension(2500) :: pendol_e_euler_t, pendol_e_euler_angle, pendol_e_euler_velocitat, &
      pendol_e_adam_t, pendol_e_adam_angle, pendol_e_adam_velocitat
   double precision, dimension(2100) :: pendol_pos_t, pendol_pos_angle, pendol_pos_velocitat, &
      pendol_neg_t, pendol_neg_angle, pendol_neg_velocitat
   double precision, dimension(300) :: pendol_300_t, pendol_300_angle, pendol_300_velocitat
   double precision, dimension(1000) :: pendol_1000_t, pendol_1000_angle, pendol_1000_velocitat
   double precision, dimension(2200) :: pendol_2200_t, pendol_2200_angle, pendol_2200_velocitat
   double precision, dimension(14500) :: pendol_14500_t, pendol_14500_angle, pendol_14500_velocitat


   common /pendol_consts/G, L, m
   double precision :: G = 3.71d0, L = 0.45d0, m = 0.510d0
   common /consts/PI
   double precision :: PI = 4.d0 * datan(1.d0)
   integer :: i

   w_N = dsqrt(G/L)
   T_N = 2.d0*PI /w_N

   open(1, file="P7-23-24-res.dat")

   ! Apartat a)
   write(1, *) "# euler pendol simplificat"
   write(1, *) "# t, angle, v_angular"
   call euler2(pendol_simplificat, 0.d0, 6*T_N, 0.0d0, 0.02d0, 1300, &
      pendol_s_t, pendol_s_velocitat, pendol_s_angle)
   do i = 1, 1300
      write(1, "(E20.12, E20.12, E20.12)") pendol_s_t(i), pendol_s_angle(i), pendol_s_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# adam pendol simplificat"
   write(1, *) "# t, angle, v_angular"
   call adams_bashforth2(pendol_simplificat, 0.d0, 6*T_N, &
      0.0d0, 0.02d0, 1300, pendol_s_t, pendol_s_velocitat, pendol_s_angle)
   do i = 1, 1300
      write(1, "(E20.12, E20.12, E20.12)") pendol_s_t(i), pendol_s_angle(i), pendol_s_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat b)
   write(1, *) "# euler pendol simple"
   write(1, *) "# t, angle, v_angular"
   call euler2(pendol, 0.d0, 6*T_N, 0.0d0, PI - 0.025d0, 1800, &
      pendol_t, pendol_velocitat, pendol_angle)
   do i = 1, 1800
      write(1, "(E20.12, E20.12, E20.12)") pendol_t(i), pendol_angle(i), pendol_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# adam pendol simple"
   write(1, *) "# t, angle, v_angular"
   call adams_bashforth2(pendol, 0.d0, 6*T_N, &
      0.0d0, PI - 0.025d0, 1800, pendol_t, pendol_velocitat, pendol_angle)
   do i = 1, 1800
      write(1, "(E20.12, E20.12, E20.12)") pendol_t(i), pendol_angle(i), pendol_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat c)
   write(1, *) "# energies 1"
   write(1, *) "# t, K, V"
   call euler2(pendol, 0.d0, 6*T_N, 0d0, 1d0, 2500, &
      pendol_e_euler_t, pendol_e_euler_velocitat, pendol_e_euler_angle)
   call adams_bashforth2(pendol, 0.d0, 6*T_N, 0d0, 1d0, 2500, &
      pendol_e_adam_t, pendol_e_adam_velocitat, pendol_e_adam_angle)
   do i = 1, 2500
      write(1, "(E20.12, E20.12, E20.12, E20.12, E20.12)") &
         pendol_e_euler_t(i), Ecine(pendol_e_euler_velocitat(i)), Epoten(pendol_e_euler_angle(i)), &
         Ecine(pendol_e_adam_velocitat(i)), Epoten(pendol_e_adam_angle(i))
   end do
   write(1, *)
   write(1, *)

   call euler2(pendol, 0.d0, 6*T_N, 0d0, PI - 0.042d0, 2500, &
      pendol_e_euler_t, pendol_e_euler_velocitat, pendol_e_euler_angle)
   call adams_bashforth2(pendol, 0.d0, 6*T_N, 0d0, PI - 0.042d0, 2500, &
      pendol_e_adam_t, pendol_e_adam_velocitat, pendol_e_adam_angle)
   write(1, *) "# energies pi"
   write(1, *) "# t, K_euler, V_euler, K_adams, V_adams"
   do i = 1, 2500
      write(1, "(E20.12, E20.12, E20.12, E20.12, E20.12)") &
         pendol_e_euler_t(i), Ecine(pendol_e_euler_velocitat(i)), Epoten(pendol_e_euler_angle(i)), &
         Ecine(pendol_e_adam_velocitat(i)), Epoten(pendol_e_adam_angle(i))
   end do
   write(1, *)
   write(1, *)

   ! Apartat d)
   call adams_bashforth2(pendol, 0d0, 7*T_N, 2*dsqrt(G/L)+0.04d0, 0d0, 2100, &
      pendol_pos_t, pendol_pos_velocitat, pendol_pos_angle)
   call adams_bashforth2(pendol, 0d0, 7*T_N, 2*dsqrt(G/L)-0.04d0, 0d0, 2100, &
      pendol_neg_t, pendol_neg_velocitat, pendol_neg_angle)
   write(1, *) "# energies pi"
   write(1, *) "# t, ang_+, v_ang_+, ang_-, v_ang_-"
   do i = 1, 2100
      write(1, "(E20.12, E20.12, E20.12, E20.12, E20.12, E20.12)") &
         pendol_pos_t(i), pendol_pos_angle(i), pendol_pos_velocitat(i), &
         pendol_neg_t(i), pendol_neg_angle(i), pendol_neg_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat e)
   write(1, *) "# convergencia 300"
   write(1, *) "# t, ang, v_ang"
   call adams_bashforth2(pendol, 0d0, 12*T_N, 0.1d0, 2.1d0, 300, &
      pendol_300_t, pendol_300_velocitat, pendol_300_angle)
   do i = 1, 300
      write(1, "(E20.12, E20.12, E20.12)") &
         pendol_300_t(i), pendol_300_angle(i), pendol_300_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 1000"
   write(1, *) "# t, ang, v_ang"
   call adams_bashforth2(pendol, 0d0, 12*T_N, 0.1d0, 2.1d0, 1000, &
      pendol_1000_t, pendol_1000_velocitat, pendol_1000_angle)
   do i = 1, 1000
      write(1, "(E20.12, E20.12, E20.12)") &
         pendol_1000_t(i), pendol_1000_angle(i), pendol_1000_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 2200"
   write(1, *) "# t, ang, v_ang"
   call adams_bashforth2(pendol, 0d0, 12*T_N, 0.1d0, 2.1d0, 2200, &
      pendol_2200_t, pendol_2200_velocitat, pendol_2200_angle)
   do i = 1, 2200
      write(1, "(E20.12, E20.12, E20.12)") &
         pendol_2200_t(i), pendol_2200_angle(i), pendol_2200_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 14500s"
   write(1, *) "# t, ang, v_ang"
   call adams_bashforth2(pendol, 0d0, 12*T_N, 0.1d0, 2.1d0, 14500, &
      pendol_14500_t, pendol_14500_velocitat, pendol_14500_angle)
   do i = 1, 14500
      write(1, "(E20.12, E20.12, E20.12)") &
         pendol_14500_t(i), pendol_14500_angle(i), pendol_14500_velocitat(i)
   end do
   write(1, *)
   write(1, *)




end program name

double precision function pendol(x, y, dy) result(retval)
   implicit none
   double precision, intent(in) :: x, y, dy

   common /pendol_consts/G, L, m
   double precision :: G, L, m
   common /consts/PI
   double precision :: PI

   retval = -G*dsin(y)/L
end function pendol

double precision function pendol_simplificat(x, y, dy) result(retval)
   implicit none
   double precision, intent(in) :: x, y, dy

   common /pendol_consts/G, L, m
   double precision :: G, L, m

   retval = -G*y/L
end function pendol_simplificat

double precision function Ecine(v_ang) result(retval)
   implicit none
   double precision, intent(in) :: v_ang
   common /pendol_consts/G, L, m
   double precision :: G, L, m
   retval = 1.d0/2.d0*m*v_ang**2*L**2
end function Ecine

double precision function Epoten(ang) result(retval)
   implicit none
   double precision, intent(in) :: ang
   common /pendol_consts/G, L, m
   double precision :: G, L, m
   retval = -m*G*L*dcos(ang)
end function Epoten

subroutine euler1(derivative, xini, xend, yini, steps, result_xs, result_ys)
   implicit none

   double precision, intent(in) :: xini, xend, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_ys

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

   result_xs(1) = xini + h
   result_ys(1) = current_y + h* derivative(current_x, current_y)

   do i = 2, steps
      result_xs(i) = xini + h*i
      result_ys(i) = current_y + 2*h* derivative(current_x, current_y)
      current_x = result_xs(i-1)
      current_y = result_ys(i-1)

      ! Codigo sucio
      write(*, *) i
   end do

end subroutine euler1

subroutine euler_doble_derivada(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
   implicit none

   double precision, intent(in) :: xini, xend, dyini, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ys

   interface
      double precision function double_derivative(x, y, dy)
         implicit none
         double precision, intent(in) :: x, y, dy
      end function double_derivative
   end interface

   integer :: i
   double precision :: h, current_x, current_y, current_dy

   current_x = xini
   current_y = yini
   current_dy = dyini

   h = (xend - xini) / dble(steps)

   do i = 1, steps
      result_xs(i) = xini + h*i
      result_dys(i) = current_dy + h*double_derivative(current_x, current_y, current_dy)
      result_ys(i) = current_y + h* current_dy
      current_x = result_xs(i)
      current_y = result_ys(i)
      current_dy = result_dys(i)

      ! Codigo sucio
      write(*, *) i
   end do
end subroutine euler_doble_derivada

subroutine adams_bashforth1(derivative, xini, xend, yini, steps, result_xs, result_ys)
   implicit none

   double precision, intent(in) :: xini, xend, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_ys

   interface
      double precision function derivative(x, y)
         implicit none
         double precision, intent(in) :: x, y
      end function derivative
   end interface

   integer :: i
   double precision :: h

   h = (xend - xini) / dble(steps)

   result_xs(1) = xini + h
   result_ys(1) = yini + h* derivative(xini, yini)

   result_xs(2) = xini + 2*h
   result_ys(2) = result_ys(1) &
      - 1.d0/2.d0*h*derivative(xini, yini) &
      + 3.d0/2.d0 * h*derivative(result_xs(1), result_ys(1))

   do i = 3, steps
      result_xs(i) = xini + h*i
      result_ys(i) = result_ys(i-1) &
         - 1.d0/2.d0*h*derivative(result_xs(i-2), result_ys(i-2)) &
         + 3.d0/2.d0*h*derivative(result_xs(i-1), result_ys(i-1))

      ! Codigo sucio
      write(*, *) i
   end do

end subroutine adams_bashforth1

subroutine adams_bashforth2(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
   implicit none
   ! todo: we may increase performance by caching the values of double derivative, in order to not compute them twice...
   double precision, intent(in) :: xini, xend, dyini, yini
   integer, intent(in) :: steps
   double precision, dimension(steps), intent(out) :: result_xs, result_dys, result_ys

   interface
      double precision function double_derivative(x, y, dy)
         implicit none
         double precision, intent(in) :: x, y, dy
      end function double_derivative
   end interface

   integer :: i
   double precision :: h

   h = (xend - xini) / dble(steps)

   ! first step using euler
   result_xs(1) = xini + h
   result_dys(1) = dyini + h*double_derivative(xini, yini, dyini)
   result_ys(1) = yini + h*result_dys(1)

   ! precalculating second step (i don't store xini, yini,
   ! dyini, therefore can't look for them in the arrays)
   result_xs(2) = xini + h*2
   result_dys(2) = result_dys(1) &
      - 1.d0/2.d0*h*double_derivative(xini, yini, dyini) &
      + 3.d0/2.d0*h*double_derivative(result_xs(1), result_ys(1), result_dys(1))
   result_ys(2) =  result_ys(1) &
      - 1.d0/2.d0*h*dyini &
      + 3.d0/2.d0*h*result_dys(1)

   do i = 3, steps
      result_xs(i) = xini + h*i
      result_dys(i) = result_dys(i-1) &
         - 1.d0/2.d0*h*double_derivative(result_xs(i-2), result_ys(i-2), result_dys(i-2)) &
         + 3.d0/2.d0*h*double_derivative(result_xs(i-1), result_ys(i-1), result_dys(i-1))
      result_ys(i) = result_ys(i-1) &
         - 1.d0/2.d0*h*result_dys(i-2) &
         + 3.d0/2.d0*h*result_dys(i-1)

      ! Codigo sucio
      write(*, *) i
   end do

end subroutine adams_bashforth2
