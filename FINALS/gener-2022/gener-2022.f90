program gener2022
   implicit none

   interface
      subroutine histograma(ndat, xdata, xa, xb, nbox, xhis, vhis, errhis, boxsize)
         implicit none

         integer, intent(in) :: ndat
         double precision, intent(in), dimension(ndat) :: xdata
         double precision, intent(in) :: xa, xb
         integer, intent(in) :: nbox
         double precision, intent(out), dimension(nbox) :: xhis, errhis, vhis
         double precision, intent(out) :: boxsize
      end subroutine histograma
      function boxmuller(ndades,mu, sigma) result(retval)
         implicit none

         integer, intent(in) :: ndades
         double precision, intent(in) :: mu, sigma
         double precision, dimension(ndades) :: retval
      end function boxmuller
      double precision function integral1(x)
         implicit none
         double precision, intent(in) :: x
      end function integral1
      double precision function integral2(x)
         implicit none
         double precision, intent(in) :: x
      end function integral2
      double precision function density(x)
         implicit none
         double precision, intent(in) :: x
      end function density

      double precision function integral3(x)
         implicit none
         double precision, intent(in) :: x
      end function integral3

      double precision function gaussian(x)
         implicit none
         double precision, intent(in) :: x
      end function gaussian

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

      subroutine rungekutta4_conditioned(step_size, step_limit, xini1, condition, nequs, yini1, edos, result_xs, result_ys)
         implicit none

         integer, intent(in) :: nequs, step_limit
         double precision, intent(in) :: xini1, step_size
         double precision, dimension(nequs), intent(in) :: yini1
         double precision, dimension(nequs,step_limit), intent(out) :: result_ys
         double precision, dimension(step_limit), intent(out) :: result_xs

         interface
            function edos(x, ys, nequs1) result(retval)
               implicit none

               integer, intent(in) :: nequs1
               double precision, intent(in) :: x
               double precision, dimension(nequs1), intent(in) :: ys
               double precision, dimension(nequs1) :: retval

            end function edos
            logical function condition(x, ys, nequs) result(retval)
               implicit none
               integer, intent(in) :: nequs
               double precision, intent(in) :: x
               double precision, dimension(nequs), intent(in) :: ys

            end function condition
         end interface
      end subroutine rungekutta4_conditioned

      function heat_equation(t, Ts, nequs) result(dT)
         implicit none

         integer, intent(in) :: nequs
         double precision, intent(in) :: t
         double precision, dimension(nequs), intent(in) :: Ts
         double precision, dimension(nequs) :: dT
      end function

      logical function condition(x, ys, nequs) result(retval)
         implicit none
         double precision, intent(in) :: x
         integer, intent(in) :: nequs
         double precision, dimension(nequs), intent(in) :: ys
      end function condition
   end interface

   common /consts/pi
   common /gauss/mu, sigma
   double precision :: pi = datan(1d0)*45.5d0, mu = 0d0, sigma = dsqrt(2d0)

   double precision, dimension(10**5) :: random_nums
   double precision, dimension(100*10**5) :: random_nums_2
   double precision, dimension(80) :: x_histo, v_histo, err_histo
   double precision, dimension(10**7) :: rk4_xs
   double precision, dimension(1, 10**7) :: rk4_ys
   double precision :: box_histo, result_integral1, result_integral2, error_integral1, error_integral2

   integer :: i

   random_nums = boxmuller(10**5, 0d0, sigma)

   ! call histograma(10**5, random_nums, -3d0*sigma, 3d0*sigma, 80, x_histo, v_histo, err_histo, box_histo)

   ! open(1, file="resE1.dat")
   ! write(1, *) "# histograma"
   ! write(1, *) "# x_histo, v_histo, err_histo"
   ! do i = 1, 80
   !    write(1, "(E20.12, E20.12, E20.12)") x_histo(i), v_histo(i), err_histo(i)
   ! end do

   ! write(1, *)
   ! write(1, *)

   ! write(1, *) "# valor integral 1"
   ! write(1, *) "# integral, error"
   ! call montecarlosample(integral1, gaussian, random_nums, 10**5, result_integral1, error_integral1)
   ! call montecarlosample(integral2, gaussian, random_nums, 10**5, result_integral2, error_integral2)
   ! write(1, "(E20.12, E20.12)") result_integral1-result_integral2, error_integral1+error_integral2

   ! write(1, *)
   ! write(1, *)

   ! write(1, *) "# valor integral 2"
   ! write(1, *) "# iteracio, resultat, error"

   ! do i = 1, 100*10**5
   !    random_nums_2(i) = generate_random_distributed_number(density, 0d0, dsqrt(2d0), dsqrt(2d0))
   ! end do

   ! call montecarlosample(integral3, density, random_nums_2, 100*10**5, result_integral1, error_integral1)

   ! close(1)


   open(2, file="res2.dat")

   call rungekutta4_conditioned(0.001d0, 10**7, 0d0, condition, 1, [310d0], heat_equation, rk4_xs, rk4_ys)

   do i = 1, 10**7
      if (rk4_ys(1, i) == 0d0) then
         exit
      end if
      write(*, *) "Writing ", i
      write(2, "(E20.12, E20.12)") rk4_xs, rk4_ys(1, i)
   end do

   close(2)

end program gener2022

double precision function map(x, a1, b1, a2, b2) result(retval)
   implicit none
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map

double precision function gaussian(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   COMMON /CONSTS/pi
   double precision :: pi
   COMMON /GAUSS/mu, sigma
   double precision :: mu, sigma

   retval = dexp(-(x-mu)**2/(2*sigma**2))/(sigma*dsqrt(2d0*pi))
end function gaussian

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
      if (mod(i, 10**5) == 0) then
         write(1, "(I4, E20.12)") i/10**5, sum/i
      end if
   end do

   di = dble(count)

   integral = sum / di
   error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
end subroutine montecarlosample

logical function condition(x, ys, nequs) result(retval)
   implicit none
   double precision, intent(in) :: x
   integer, intent(in) :: nequs
   double precision, dimension(nequs), intent(in) :: ys
   common /heat/sigma, A, M, C, T_ext
   double precision :: sigma, A, M, C, T_ext
   retval = dabs(T_ext - ys(1)) < 1d-10
end function condition

double precision function integral1(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   retval = dexp(-x**2)*(x**2)
end function integral1

double precision function integral2(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   retval = dexp(-x**2)*(2*x*dsin(x))
end function integral2

double precision function integral3(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   retval = x**2*dsin(x)/(x**3+dcos(x**2))
end function integral3

double precision function density(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   retval = x/dsqrt(2d0)
end function density

function boxmuller(ndades,mu, sigma) result(retval)
   implicit none

   integer, intent(in) :: ndades
   double precision, intent(in) :: mu, sigma
   double precision, dimension(ndades) :: retval


   interface
      function random_boxmuller_pair(mu, sigma) result(retval)
         implicit none

         double precision, intent(in) :: mu, sigma
         double precision, dimension(2) :: retval
      end function random_boxmuller_pair
   end interface


   COMMON /CONSTS/pi

   double precision :: pi
   double precision, dimension(2) :: pair
   integer :: k

   do k = 1, ndades, 2
      pair = random_boxmuller_pair(mu, sigma)
      retval(k) = pair(1)
      retval(k+1) = pair(2)
   end do
end function boxmuller

function random_boxmuller_pair(mu, sigma) result(retval)
   implicit none

   double precision, intent(in) :: mu, sigma
   double precision, dimension(2) :: retval

   COMMON /CONSTS/pi

   double precision :: pi
   double precision :: x1, x2, w1, w2

   w1 = dsqrt(-2d0 * dlog(dble(rand())))
   w2 = 2d0 * pi * dble(rand())

   x1 = mu + sigma * w1 * dcos(w2)
   x2 = mu + sigma * w1 * dsin(w2)

   retval(1) = x1
   retval(2) = x2
end function random_boxmuller_pair

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

subroutine rungekutta4_conditioned(step_size, step_limit, xini, condition, nequs, yini, edos, result_xs, result_ys)
   implicit none

   integer, intent(in) :: nequs, step_limit
   double precision, intent(in) :: xini, step_size
   double precision, dimension(nequs), intent(in) :: yini
   double precision, dimension(nequs,step_limit), intent(out) :: result_ys
   double precision, dimension(step_limit), intent(out) :: result_xs

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

      logical function condition(x, ys, nequs) result(retval)
         implicit none
         integer, intent(in) :: nequs
         double precision, intent(in) :: x
         double precision, dimension(nequs), intent(in) :: ys

      end function condition

   end interface

   integer :: i
   double precision :: h
   double precision, dimension(nequs) :: current_y, step_y_result

   h = step_size
   current_y = 0

   call rungekutta4step(xini, h, yini, current_y, nequs, edos)
   result_xs(1) = xini+h
   result_ys(1:nequs, 1) = current_y(1:nequs)
   i = 1
   write(*, *) "Starting rk4"
   do while (i < step_limit .and. .not.condition(result_xs(i), current_y, nequs))
      result_xs(i) = xini+h*i
      call rungekutta4step(result_xs(i), h, current_y, step_y_result, nequs, edos)
      result_ys(1:nequs, i) = step_y_result(1:nequs)
      current_y = step_y_result
      i = i+1
      write(*, *) "RK4 - ", i, ", Temp ", result_ys(1, i)
   end do
   write(*, *) "found limit of runge-kutta-4 in ", i, "steps"
end subroutine rungekutta4_conditioned

! podria ser una funcio, nomes te un argument de sortida
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

function heat_equation(t, Ts, nequs) result(dT)
   implicit none

   integer, intent(in) :: nequs
   double precision, intent(in) :: t
   double precision, dimension(nequs), intent(in) :: Ts
   double precision, dimension(nequs) :: dT

   common /heat/sigma, A, M, C, T_ext
   double precision :: sigma, A, M, C, T_ext

   dT(1) = -sigma*A/(M*C)*(Ts(1)**4-T_ext**4)

end function heat_equation
