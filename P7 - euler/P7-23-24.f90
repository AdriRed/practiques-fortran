program practica7
   implicit none

   interface
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
      end subroutine euler_doble_derivada
      subroutine euler2_doble_derivada(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
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
      end subroutine euler2_doble_derivada

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
   double precision, dimension(1500) :: pendol_s_t, pendol_s_angle, pendol_s_velocitat
   double precision, dimension(1500) :: pendol_t, pendol_angle, pendol_velocitat
   double precision, dimension(1500) :: pendol_e_euler_t, pendol_e_euler_angle, pendol_e_euler_velocitat, &
      pendol_e_euler2_t, pendol_e_euler2_angle, pendol_e_euler2_velocitat
   double precision, dimension(6000) :: pendol_pos_t, pendol_pos_angle, pendol_pos_velocitat, &
      pendol_neg_t, pendol_neg_angle, pendol_neg_velocitat
   double precision, dimension(300) :: pendol_300_t, pendol_300_angle, pendol_300_velocitat
   double precision, dimension(550) :: pendol_550_t, pendol_550_angle, pendol_550_velocitat
   double precision, dimension(1000) :: pendol_1000_t, pendol_1000_angle, pendol_1000_velocitat
   double precision, dimension(20000) :: pendol_20000_t, pendol_20000_angle, pendol_20000_velocitat


   common /pendol_consts/G, L, m
   double precision :: G = 10.44d0, L = 1.07d0, m = 0.98d0
   common /consts/PI
   double precision :: PI = 4.d0 * datan(1.d0)
   integer :: i

   w_N = dsqrt(G/L)
   T_N = 2.d0*PI /w_N

   open(1, file="P7-23-24-res.dat")

   ! Apartat a)
   write(1, *) "# euler pendol simplificat"
   write(1, *) "# t, angle, v_angular"
   call euler_doble_derivada(pendol_simplificat, 0.d0, 7*T_N, 0.0d0, 0.025d0, 1500, &
      pendol_s_t, pendol_s_velocitat, pendol_s_angle)
   do i = 1, 1500
      write(1, "(E20.12, E20.12, E20.12)") pendol_s_t(i), pendol_s_angle(i), pendol_s_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# euler2 pendol simplificat"
   write(1, *) "# t, angle, v_angular"
   call euler2_doble_derivada(pendol_simplificat, 0.d0, 7*T_N, &
      0.0d0, 0.025d0, 1500, pendol_s_t, pendol_s_velocitat, pendol_s_angle)
   do i = 1, 1500
      write(1, "(E20.12, E20.12, E20.12)") pendol_s_t(i), pendol_s_angle(i), pendol_s_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat b)
   write(1, *) "# euler pendol simple"
   write(1, *) "# t, angle, v_angular"
   call euler_doble_derivada(pendol, 0.d0, 7*T_N, 0.0d0, PI - 0.15d0, 1500, &
      pendol_t, pendol_velocitat, pendol_angle)
   do i = 1, 1500
      write(1, "(E20.12, E20.12, E20.12)") pendol_t(i), pendol_angle(i), pendol_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# euler2 pendol simple"
   write(1, *) "# t, angle, v_angular"
   call euler2_doble_derivada(pendol, 0.d0, 7*T_N, &
      0.0d0, PI - 0.15d0, 1500, pendol_t, pendol_velocitat, pendol_angle)
   do i = 1, 1500
      write(1, "(E20.12, E20.12, E20.12)") pendol_t(i), pendol_angle(i), pendol_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat c)
   call euler_doble_derivada(pendol, 0.d0, 7*T_N, 0d0, PI - 0.025d0, 1500, &
      pendol_e_euler_t, pendol_e_euler_velocitat, pendol_e_euler_angle)
   call euler2_doble_derivada(pendol, 0.d0, 7*T_N, 0d0, PI - 0.025d0, 1500, &
      pendol_e_euler2_t, pendol_e_euler2_velocitat, pendol_e_euler2_angle)
   write(1, *) "# energies"
   write(1, *) "# t, K_euler, V_euler, K_euler2, V_euler2"
   do i = 1, 1500
      write(1, "(E20.12, E20.12, E20.12, E20.12, E20.12)") &
         pendol_e_euler_t(i), Ecine(pendol_e_euler_velocitat(i)), Epoten(pendol_e_euler_angle(i)), &
         Ecine(pendol_e_euler2_velocitat(i)), Epoten(pendol_e_euler2_angle(i))
   end do
   write(1, *)
   write(1, *)

   ! Apartat d)
   call euler2_doble_derivada(pendol, 0d0, 7*T_N, 2*dsqrt(G/L)+0.05d0, 0d0, 6000, &
      pendol_pos_t, pendol_pos_velocitat, pendol_pos_angle)
   call euler2_doble_derivada(pendol, 0d0, 7*T_N, 2*dsqrt(G/L)-0.05d0, 0d0, 6000, &
      pendol_neg_t, pendol_neg_velocitat, pendol_neg_angle)
   write(1, *) "# transicio"
   write(1, *) "# t, ang_+, v_ang_+, ang_-, v_ang_-"
   do i = 1, 6000
      write(1, "(E20.12, E20.12, E20.12, E20.12, E20.12, E20.12)") &
         pendol_pos_t(i), pendol_pos_angle(i), pendol_pos_velocitat(i), &
         pendol_neg_t(i), pendol_neg_angle(i), pendol_neg_velocitat(i)
   end do
   write(1, *)
   write(1, *)

   ! Apartat e)
   write(1, *) "# convergencia 300"
   write(1, *) "# t, ang, v_ang"
   call euler2_doble_derivada(pendol, 0d0, 11*T_N, 0.d0, 2.87d0, 300, &
      pendol_300_t, pendol_300_velocitat, pendol_300_angle)
   do i = 1, 300
      write(1, "(E20.12, E20.12)") &
         pendol_300_t(i), Epoten(pendol_300_angle(i)) + Ecine(pendol_300_velocitat(i))
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 550"
   write(1, *) "# t, ang, v_ang"
   call euler2_doble_derivada(pendol, 0d0, 11*T_N, 0.d0, 2.87d0, 550, &
      pendol_550_t, pendol_550_velocitat, pendol_550_angle)
   do i = 1, 550
      write(1, "(E20.12, E20.12)") &
         pendol_550_t(i), Epoten(pendol_550_angle(i)) + Ecine(pendol_550_velocitat(i))
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 1000"
   write(1, *) "# t, ang, v_ang"
   call euler2_doble_derivada(pendol, 0d0, 11*T_N, 0.d0, 2.87d0, 1000, &
      pendol_1000_t, pendol_1000_velocitat, pendol_1000_angle)
   do i = 1, 1000
      write(1, "(E20.12, E20.12)") &
         pendol_1000_t(i), Epoten(pendol_1000_angle(i)) + Ecine(pendol_1000_velocitat(i))
   end do
   write(1, *)
   write(1, *)

   write(1, *) "# convergencia 20000"
   write(1, *) "# t, ang, v_ang"
   call euler2_doble_derivada(pendol, 0d0, 11*T_N, 0.d0, 2.87d0, 20000, &
      pendol_20000_t, pendol_20000_velocitat, pendol_20000_angle)
   do i = 1, 20000
      write(1, "(E20.12, E20.12)") &
         pendol_20000_t(i), Epoten(pendol_20000_angle(i)) + Ecine(pendol_20000_velocitat(i))
   end do
   write(1, *)
   write(1, *)


   ! Apartat extra)
   write(1, *) "# animacio"
   write(1, *) "# t, angle, v_angular"
   call euler2_doble_derivada(pendol, 0.d0, 7*T_N, &
      0.0d0, PI - 0.15d0, 1000, pendol_1000_t, pendol_1000_velocitat, pendol_1000_angle)
   do i = 1, 1000
      write(1, "(E20.12, E20.12, E20.12)") pendol_1000_t(i), pendol_1000_angle(i), pendol_1000_velocitat(i)
      write(1, *)
      write(1, *)   
   end do
   write(1, *)
   write(1, *)

end program practica7

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
   end do
end subroutine euler_doble_derivada

subroutine euler2_doble_derivada(double_derivative, xini, xend, dyini, yini, steps, result_xs, result_dys, result_ys)
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
   double precision :: h, current_x, current_y, current_dy, k1, k2

   current_x = xini
   current_y = yini
   current_dy = dyini

   h = (xend - xini) / dble(steps)

   ! euler first order way
   result_xs(1) = xini + h
   result_dys(1) = dyini + h*double_derivative(xini, yini, dyini)
   result_ys(1) = yini + h*result_dys(1)

   do i = 2, steps
      result_xs(i) = xini + h*i
      k1 = double_derivative(result_xs(i-1), result_ys(i-1), result_dys(i-1))                         !--------------------!
      k2 = double_derivative(result_xs(i-1) + 3.d0*h/4.d0, result_ys(i-1) + 3.d0*h/4.d0*k1, result_dys(i-1)+ 3.d0*h/4.d0*k1)
      result_dys(i) = result_dys(i-1) &
         + 1.d0/3.d0*h*k1 &
         + 2.d0/3.d0*h*k2
      result_ys(i) = result_ys(i-1) &
         + h*result_dys(i)
   end do
end subroutine euler2_doble_derivada
