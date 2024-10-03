program prepra3
   implicit none

   interface
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


      end subroutine newtonraphson
      subroutine sinh_polinomi(x, f_x, df_x)      
         double precision, intent(in) :: x
         double precision, intent(out) :: f_x, df_x
      end subroutine sinh_polinomi
      subroutine derivataula(length, xvalues, f_xvalues, df_xvalues)
         integer, intent(in) :: length
         double precision, dimension(length), intent(in) :: xvalues, f_xvalues
         double precision, dimension(length), intent(out) :: df_xvalues
      end subroutine derivataula
   end interface

   common /CONSTS/PI

   double precision :: PI = 3.14159265d0
   double precision :: E, step = 0.0001, f_e, df_e, arrel
   double precision, dimension(25) :: Es25, f_Es25, df_Es25_Approx, df_Es25
   double precision, dimension(230) :: Es230, f_Es230, df_Es230_Approx, df_Es230
   double precision, dimension(9) :: E_0
   integer :: iteracions, i 
   iteracions = 10
   
   !Apartat 2.a
   E = 0
   open(1, file='P3-23-24-res.dat')
   write(1,  *) "# Dades de F(E) i dF(E)"

   do while (E <= 2*PI)
      call sinh_polinomi(E, f_e, df_e)
      write(1, "(E20.12, E20.12, E20.1 2)") E, f_e, df_e
      E = E+step
   end do

   !Apartat 2.b
   write(*, *) "Métode de bisecció"
   call biseccio(sinh_polinomi, 0.5d0, 0.8d0, 1.d-12, iteracions, arrel)
   write(*, "(a, E20.12, a, I3, a)") "Arrel de la funció del P(E) trobada a E = ", arrel, " amb ", iteracions, " iteracions" 

   call biseccio(sinh_polinomi, 0.8d0, 1.4d0, 1.d-12, iteracions, arrel)
   write(*, "(a, E20.12, a, I3, a)") "Arrel de la funció del P(E) trobada a E = ", arrel, " amb ", iteracions, " iteracions" 

   call biseccio(sinh_polinomi, 1.4d0, 1.5d0, 1.d-12, iteracions, arrel)
   write(*, "(a, E20.12, a, I3, a)") "Arrel de la funció del P(E) trobada a E = ", arrel, " amb ", iteracions, " iteracions" 

   !Apartat 2.c
   write(1, *)
   write(1, *)
   write(1,  *) "# Càlcul d'arrels mitjançant Newton-Raphson"

   write(*, *)

   E_0 = (/ 0.1, 0.2, 0.65, 0.7, 1.3, 2.4, 2.6, 3.9, 5.3 /)
   do i = 1, 9
      E = E_0(i)
      call newtonraphson(sinh_polinomi, E, 1.d-12, iteracions, arrel)
      write(*, "(a, E20.12, a, E20.12, a, I3, a)") &
         "Arrel de la funcio del P(E) trobada a E = ", arrel, " des de E_0 = ", E, " amb ", iteracions, " iteracions" 
      write(1, "(E20.12, I3)") E, iteracions
   end do
   close(1)

   !Apartat 3
   open(2, file="P3-23-24-res3-n25.dat")
   step = 2d0/25d0*PI
   E = 0d0
   do i = 1, 25
      call sinh_polinomi(E, f_e, df_e)
      Es25(i) = E
      f_Es25(i) = f_e
      df_Es25(i) = df_e
      E = E+step
   end do
   call derivataula(25, Es25, f_Es25, df_Es25_Approx)
   do i = 1, 25
      write(2, "(E20.12, E20.12, E20.12, E20.12, E20.12)") Es25(i), f_Es25(i), df_Es25_Approx(i), df_Es25(i)
   end do
   close(2)

   open(3, file="P3-23-24-res3-n230.dat")
   step = 2d0/230d0*PI
   E = 0d0
   do i = 1, 230
      call sinh_polinomi(E, f_e, df_e)
      Es230(i) = E
      f_Es230(i) = f_e
      df_Es230(i) = df_e
      E = E+step
   end do
   call derivataula(230, Es230, f_Es230, df_Es230_Approx)
   do i = 1, 230
      write(3, "(E20.12, E20.12, E20.12, E20.12, E20.12)") Es230(i), f_Es230(i), df_Es230_Approx(i), df_Es230(i)
   end do
   close(3)




end program prepra3

subroutine sinh_polinomi(x, f_x, df_x)
   interface
      subroutine sinh_and_derivative(x1, f_x1, df_x1)
         double precision, intent(in) :: x1
         double precision, intent(out) :: f_x1, df_x1
      end subroutine sinh_and_derivative
      subroutine test_polinomi(x2, f_x2, df_x2)
         double precision, intent(in) :: x2
         double precision, intent(out) :: f_x2, df_x2
      end subroutine test_polinomi
      subroutine mult_funcs(x12, f, g, fg_x, dfg_x)
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

         double precision, intent(in) :: x12
         double precision, intent(out) :: fg_x, dfg_x
      end subroutine mult_funcs

   end interface

   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   call mult_funcs(x, sinh_and_derivative, test_polinomi, f_x, df_x)

end subroutine sinh_polinomi

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

subroutine sinh_and_derivative(x, f_x, df_x)
   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x

   f_x = sinh(x)
   df_x = cosh(x)

end subroutine sinh_and_derivative

subroutine test_polinomi(x, f_x, df_x)

   COMMON /CONSTS/PI

   double precision, intent(in) :: x
   double precision, intent(out) :: f_x, df_x
   double precision :: PI

   f_x = -57./160. * PI + (57./80. + 17./20. * PI)*x - (17./10. + PI/2.)*x**2 + x**3
   df_x = (57./80. + 17./20. * PI) - 2*(17./10. + PI/2.)*x + 3*x**2
end subroutine test_polinomi

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