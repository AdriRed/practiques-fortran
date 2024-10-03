! Adria Rojo Merida

program final
   implicit none

   interface
      subroutine rungekutta4(steps1, xini1, xend1, nequs, yini1, edos, result_xs, result_ys)
         implicit none

         integer, intent(in) :: nequs, steps1
         double precision, intent(in) :: xini1, xend1
         double precision, dimension(nequs), intent(in) :: yini1
         double precision, dimension(nequs,steps1), intent(out) :: result_ys
         double precision, dimension(steps1), intent(out) :: result_xs

         interface
            function edos(x, ys1, nequs1) result(retval)
               implicit none

               integer, intent(in) :: nequs1
               double precision, intent(in) :: x
               double precision, dimension(nequs1), intent(in) :: ys1
               double precision, dimension(nequs1) :: retval

            end function edos
         end interface
      end subroutine rungekutta4
      function problema1(t, ys1, nequs) result(retval)
         implicit none
         integer, intent(in) :: nequs
         double precision, intent(in) :: t
         double precision, dimension(nequs), intent(in) :: ys1
         double precision, dimension(nequs) :: retval
      end function
      double precision function energy(z, phi)
         implicit none
         double precision, intent(in) :: z, phi
      end function energy

      double precision function simpsonstresvuit(x1, x2, k1, func) result(retval)
         implicit none
         interface
            double precision function func(x)
               implicit none
               double precision, intent(in) :: x
            end function func
         end interface
         double precision, intent(in) :: x1, x2
         integer, intent(in) :: k1
      end function simpsonstresvuit
      double precision function integral1(x)
         implicit none
         double precision, intent(in) :: x
      end function integral1
      function random_boxmuller_pair(mu, sigma) result(retval)
         implicit none

         double precision, intent(in) :: mu, sigma
         double precision, dimension(2) :: retval
      end function random_boxmuller_pair
      double precision function gaussian(x) 
           implicit none
           double precision, intent(in) :: x
      end function gaussian
   end interface

   common /p1/delta, omega
   double precision :: delta = 2.5d0, omega = 1d0
   integer :: i, k, steps = 10000
   double precision, dimension(10000) :: ts
   double precision, dimension(2,10000) :: ys
   double precision, dimension(2) :: zs

   common /integrals/x_0, gamma
   double precision :: x_0 = 1d0, gamma = 0.5d0
   common /consts/pi
   double precision :: pi
   common /gauss/mu, sigma
   double precision :: mu = 0d0, sigma = 100d0 ! seguro que cojo > 99% de los puntos con esta sigma
   double precision :: i1, err
   double precision, dimension(2) :: pair
   double precision, dimension(100000) :: random_nums

   pi = 4d0*datan(1d0)

   open(1, file="Exa-jan-24-res1.dat")
   ! Problema 1
   zs = [0.3d0, 0.9d0]

   do i = 1, 2
      write(*, *) "Calculando Runge-Kutta 4 para z(0) = ", zs(i)
      call rungekutta4(steps, 0d0, 10d0, 2, [zs(i), 0d0], problema1, ts, ys)
      write(1, *) "# Valores de z y phi para z(0) = ", zs(i)
      write(1, *) "# t, z, phi"
      do k = 1, steps
         ! como es un angulo, si encuentro alguna phi < 0 le sumo 2pi, para cumplir
         ! la restricción del enunciado, de que phi € [0, 2pi)
         ! por esta razón la figura queda "partida"
         if (ys(2, k) < 0) then
            ys(2, k) = ys(2, k) + 2*pi
         end if

         write(1, "(E20.12, E20.12, E20.12)") ts(k), ys(1, k), ys(2, k)
      end do
      write(1, *)
      write(1, *)
      write(1, *) "# Energias inicial y final para z(0) = ", zs(i)
      write(1, "(E20.12, E20.12)") energy(ys(1, 1), ys(2, 1)), energy(ys(1, steps), ys(2, steps))
      write(1, *)
      write(1, *)
   end do

   ! Problema 2
   write(1, *) "# Integral tancada "
   write(1, *) "# i1"
   ! con N = 3**10 me paso, pero ya va bien
   write(1, "(E20.12)") simpsonstresvuit(0d0, pi, 10, integral1) 
   write(1, *)
   write(1, *)

   write(*, *) "Calculando valores aleatorios"
   do i = 1, 100000, 2
      pair = random_boxmuller_pair(mu, sigma)
      random_nums(i) = pair(1)
      random_nums(i+1) = pair(2)
   end do

   write(*, *) "Calculando Montcarlo"
   write(1, *) "# Integral Montecarlo"
   write(1, *) "# i, I, err"
   call montecarlosample(integral1, gaussian, random_nums, 100000, i1, err)


   close(1)

end program final

function random_boxmuller_pair(mu, sigma) result(retval)
   implicit none

   double precision, intent(in) :: mu, sigma
   double precision, dimension(2) :: retval

   COMMON /CONSTS/PI

   double precision :: PI
   double precision :: x1, x2, w1, w2

   w1 = dsqrt(-2d0 * dlog(dble(rand())))
   w2 = 2d0 * PI * dble(rand())

   x1 = mu + sigma * w1 * dcos(w2)
   x2 = mu + sigma * w1 * dsin(w2)

   retval(1) = x1
   retval(2) = x2
end function random_boxmuller_pair

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

      ! codigo sucio
      if (modulo(i, 10000) == 0) then
         di = dble(i)

         integral = sum / di
         error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
         write(*, *) "Montecarlo ", i
         write(1, "(I9, E20.12, E20.12)")  i/10000, &
            integral, error
      end if
   end do

   di = dble(count)

   integral = sum / di
   error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
end subroutine montecarlosample

double precision function gaussian(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   common /consts/pi
   double precision :: pi
   common /gauss/mu, sigma
   double precision :: mu, sigma

   retval = dexp(-(x-mu)**2/(2*sigma**2))/(sigma*dsqrt(2d0*pi))
end function gaussian

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

double precision function integral1(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /integrals/x_0, gamma
   double precision :: x_0, gamma
   retval = dsin(x)**2/((x-x_0)**2 + gamma**2)
end function integral1

double precision function energy(z, phi) result(retval)
   implicit none
   double precision, intent(in) :: z, phi
   common /p1/delta, omega
   double precision :: delta, omega

   retval = delta*z**2/2d0 - dsqrt(1-z**2)*dcos(phi)*omega
end function energy

! ys(1) = z
! ys(2) = phi
function problema1(t, ys, nequs) result(retval)
   implicit none
   integer, intent(in) :: nequs
   double precision, intent(in) :: t
   double precision, dimension(nequs), intent(in) :: ys
   double precision, dimension(nequs) :: retval

   common /p1/delta, omega
   double precision :: delta, omega


   retval(1) = -omega*dsqrt(1d0-ys(1)**2)*dsin(ys(2))
   retval(2) = delta*ys(1) + omega*ys(1)/dsqrt(1d0-ys(1)**2) * dcos(ys(2))
end function problema1


subroutine rungekutta4(steps, xini, xend, nequs, yini, edos, result_xs, result_ys)
   implicit none

   integer, intent(in) :: nequs, steps
   double precision, intent(in) :: xini, xend
   double precision, dimension(nequs), intent(in) :: yini
   double precision, dimension(nequs,steps), intent(out) :: result_ys
   double precision, dimension(steps), intent(out) :: result_xs

   interface
      function edos(x, ys, nequs1) result(retval)
         implicit none

         integer, intent(in) :: nequs1
         double precision, intent(in) :: x
         double precision, dimension(nequs1), intent(in) :: ys
         double precision, dimension(nequs1) :: retval

      end function edos
      subroutine rungekutta4step(x, h1, funcin, dfuncout, nequs1, edofuncio)
         implicit none
         integer, intent(in) :: nequs1
         double precision, dimension(nequs1), intent(in) :: funcin
         double precision, dimension(nequs1), intent(out) :: dfuncout
         double precision, intent(in) :: x, h1

         interface
            function edofuncio(x1, ys, nequs2) result(retval)
               implicit none

               integer, intent(in) :: nequs2
               double precision, intent(in) :: x1
               double precision, dimension(nequs2), intent(in) :: ys
               double precision, dimension(nequs2) :: retval

            end function edofuncio
         end interface
      end subroutine rungekutta4step
   end interface

   integer :: i
   double precision :: h
   double precision, dimension(nequs) :: current_y, step_y_result

   h = (xend - xini) / steps
   current_y = 0

   call rungekutta4step(xini, h, yini, current_y, nequs, edos)
   result_xs(1) = xini+h
   result_ys(1:nequs, 1) = current_y(1:nequs)

   do i = 2, steps
      result_xs(i) = xini+h*i
      call rungekutta4step(result_xs(i), h, current_y, step_y_result, nequs, edos)
      result_ys(1:nequs, i) = step_y_result(1:nequs)
      current_y = step_y_result
   end do
end subroutine rungekutta4

subroutine rungekutta4step(x0, h, y0, y1, nequs, edofuncio)
   implicit none
   integer, intent(in) :: nequs
   double precision, dimension(nequs), intent(in) :: y0
   double precision, dimension(nequs), intent(out) :: y1
   double precision, intent(in) :: x0, h

   interface
      function edofuncio(x, ys, nequs1) result(retval)
         implicit none

         integer, intent(in) :: nequs1
         double precision, intent(in) :: x
         double precision, dimension(nequs1), intent(in) :: ys
         double precision, dimension(nequs1) :: retval

      end function edofuncio
   end interface

   double precision, dimension(nequs) :: k1, k2, k3, k4

   k1 = edofuncio(x0, y0, nequs)
   k2 = edofuncio(x0 + h/2d0, y0 + h*k1/2d0, nequs)
   k3 = edofuncio(x0 + h/2d0, y0 + h*k2/2d0, nequs)
   k4 = edofuncio(x0 + h, y0 + h*k3, nequs)
   y1 = y0 + h/6d0 * (k1 + 2*k2 + 2*k3 + k4)

end subroutine rungekutta4step
