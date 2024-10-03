program prepra6
   implicit none

   interface
      double precision function integral_1(x) result(retval)
         double precision, intent(in) :: x
      end function integral_1
      double precision function integral_2(x) result(retval)
         double precision, intent(in) :: x
      end function integral_2
      double precision function integral_3(x) result(retval)
         double precision, intent(in) :: x
      end function integral_3
      double precision function integral_4(x) result(retval)
         double precision, intent(in) :: x
      end function integral_4
      double precision function integral_5(x) result(retval)
         double precision, intent(in) :: x
      end function integral_5
      double precision function integral_6(xs) result(retval)
         implicit none
         double precision, dimension(:), intent(in) :: xs
      end function integral_6
      double precision function func_p(x) result(retval)
         double precision, intent(in) :: x
      end function func_p
      double precision function func_g(xs) result(retval)
         implicit none
         double precision, intent(in), dimension(:) :: xs
      end function func_g
      double precision function gaussian(x) result(retval)
         double precision, intent(in) :: x
      end function gaussian

      subroutine montecarlocru(func, count, a, b, integral,  error)
         implicit none

         integer, intent(in) :: count
         double precision, intent(in) :: a, b
         double precision, intent(out) :: integral, error

         interface
            double precision function func(x1)
               double precision, intent(in) :: x1
            end function func
            double precision function map(x1, a1, b1, a2, b2) result(retval)
               double precision, intent(in) :: x1, b1, a1, a2, b2
            end function map
         end interface
      end subroutine montecarlocru

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

      subroutine montecarlosample_multidim(func, distr, dims, random_vectors, count, integral, error)
         implicit none
         integer, intent(in) :: count, dims
         double precision, intent(in), dimension(count, dims) :: random_vectors
         double precision, intent(out) :: integral, error
         interface
            double precision function func(xs)
               double precision, intent(in), dimension(:) :: xs
            end function func
            double precision function distr(xs)
               double precision, intent(in), dimension(:) :: xs
            end function distr
         end interface
      end subroutine montecarlosample_multidim

      double precision function generate_random_distributed_number(distr, a, b, cotasup) result(retval)
         double precision, intent(in) :: a, b, cotasup

         interface
            double precision function distr(x)
               double precision, intent(in) :: x
            end function distr
         end interface
      end function generate_random_distributed_number

      function random_boxmuller_pair(mu1, sigma1) result(retval)
         implicit none

         double precision, intent(in) :: mu1, sigma1
         double precision, dimension(2) :: retval
      end function

      function random_metropolis(dims, quantity, distr) result(retval)
         implicit none

         interface
            double precision function distr(xs)
               implicit none
               double precision, dimension(:), intent(in) :: xs
            end function distr
         end interface
         integer, intent(in) :: dims, quantity
         double precision, dimension(quantity,dims) :: retval
      end function random_metropolis
   end interface

   COMMON /CONSTS/PI
   double precision :: PI = datan(1d0)*4d0

   COMMON /GAUSS/MU, SIGMA
   double precision :: MU = 0d0, SIGMA = 1d0/dsqrt(2d0)

   integer :: i
   double precision, dimension(1100000) :: random_nums, box_muller_nums
   double precision, dimension(2) :: boxmuller_pair
   double precision, dimension(210000,5) :: random_vects
   double precision :: res_integral1, error1, res_integral2, error2
   double precision :: res_integral3, error3, res_integral4, error4, res_integral5, error5
   double precision :: res_integral6, error6


   call srand(20401021)

   open(1, file="data.dat")
   write(1, *) "# montecarlo cru"
   write(1, *) "# N, I_1, estimacio_1, real_1, I_2, estimacio_2, real_2"

   do i = 2000, 120000, 2000
      call montecarlocru(integral_1, i, -PI, PI, res_integral1, error1)
      call montecarlocru(integral_2, i, -2d0*PI, 2d0*PI, res_integral2, error2)

      write(1, "(I6, E20.12, E20.12, E20.12, E20.12, E20.12, E20.12)") i, &
         res_integral1, error1, &
         res_integral2, error2
   end do
   write(*, *) "MONTECARLO CRU OK"
   write(1, *)
   write(1, *)

   do i = 1, 1100000
      random_nums(i) = generate_random_distributed_number(func_p, 0d0, PI, 10d0/3d0)
   end do

   do i = 1, 1100000, 2
      boxmuller_pair = random_boxmuller_pair(MU, SIGMA)
      box_muller_nums(i) = boxmuller_pair(1)
      box_muller_nums(i+1) = boxmuller_pair(2)
   end do

   write(1, *) "# montecarlo sample I_3"
   write(1, *) "# N, I_3, err_3"

   call montecarlosample(integral_3, func_p, random_nums, 1100000, res_integral3, error3)
   write(1, *)
   write(1, *)

   write(1, *) "# montecarlo sample I_4"
   write(1, *) "# N, I_4, err_4"
   call montecarlosample(integral_4, func_p, random_nums, 1100000, res_integral4, error4)
   write(1, *)
   write(1, *)

   write(1, *) "# montecarlo sample I_5"
   write(1, *) "# N, I_5, err_5"
   call montecarlosample(integral_5, gaussian, box_muller_nums, 1100000, res_integral5, error5)
   write(1, *)
   write(1, *)

   write(*, *) "MONTECARLO SAMPLE OK"

   random_vects = random_metropolis(5, 210000, func_g)
   write(1, *) "# montecarlo multidimensional I_5"
   write(1, *) "# N, I_6, err_6"
   call montecarlosample_multidim(integral_6, func_g, 5, random_vects, 210000, res_integral6, error6)
   write(1, *)
   write(1, *)
   write(*, *) "MONTECARLO SAMPLE MULTIDIM OK"

end program prepra6

double precision function gaussian(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   COMMON /CONSTS/pi
   double precision :: pi
   COMMON /GAUSS/mu, sigma
   double precision :: mu, sigma

   retval = dexp(-(x-mu)**2/(2*sigma**2))/(sigma*dsqrt(2d0*pi))
end function gaussian

double precision function map(x, a1, b1, a2, b2) result(retval)
   implicit none
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map

subroutine montecarlocru(func, count, a, b, integral, error)
   implicit none

   integer, intent(in) :: count
   double precision, intent(in) :: a, b
   double precision, intent(out) :: integral, error

   interface
      double precision function func(x1)
         double precision, intent(in) :: x1
      end function func
      double precision function map(x1, a1, b1, a2, b2) result(retval)
         double precision, intent(in) :: x1, b1, a1, a2, b2
      end function map
   end interface

   double precision :: sum, sum2, x, dcount
   integer :: i

   sum = 0d0
   sum2 = 0d0

   do i = 1, count
      x = map(dble(rand()), 0d0, 1d0, a, b)
      sum = sum + func(x)*(b-a)
      sum2 = sum2 + (func(x)*(b-a))**2
   end do

   dcount = dble(count)

   integral = sum / dcount
   error = 1/dsqrt(dcount) * dsqrt(sum2/dcount - integral**2)
end subroutine montecarlocru

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
      if (modulo(i, 5000) == 0) then
         di = dble(i)

         integral = sum / di
         error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
         write(*, *) i
         write(1, "(I9, E20.12, E20.12)")  i, &
            integral, error
      end if
   end do

   di = dble(count)

   integral = sum / di
   error = 1/dsqrt(di) * dsqrt(sum2/di - integral**2)
end subroutine montecarlosample

subroutine montecarlosample_multidim(func, distr, dims, random_vectors, count, integral, error)
   implicit none
   integer, intent(in) :: count, dims
   double precision, intent(in), dimension(count, dims) :: random_vectors
   double precision, intent(out) :: integral, error
   interface
      double precision function func(xs)
         double precision, intent(in), dimension(:) :: xs
      end function func
      double precision function distr(xs)
         double precision, intent(in), dimension(:) :: xs
      end function distr
   end interface

   double precision :: sum, sum2, dcount
   double precision, dimension(dims) :: current_vector
   integer :: i, k

   sum = 0d0
   sum2 = 0d0

   do i = 1, count
      do k = 1, dims !random vector choosing
         current_vector(k) = random_vectors(i, k)
      end do
      sum = sum + func(current_vector)/distr(current_vector)
      sum2 = sum2 + (func(current_vector)**2/distr(current_vector)**2)

      ! codigo sucio
      if (modulo(i, 1500) == 0) then
         dcount = dble(i)

         integral = sum / dcount
         error = 1/dsqrt(dcount) * dsqrt(sum2/dcount - integral**2)
         write(*, *) i
         write(1, "(I9, E20.12, E20.12)")  i, &
            integral, error
      end if
   end do

   dcount = dble(count)

   integral = sum / dcount
   error = 1/dsqrt(dcount) * dsqrt(sum2/dcount - integral**2)
end subroutine montecarlosample_multidim

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

double precision function generate_random_distributed_number(distr, a, b, cotasup) result(retval)
   implicit none

   double precision, intent(in) :: a, b, cotasup

   interface
      double precision function distr(x1)
         double precision, intent(in) :: x1
      end function distr
      double precision function map(x1, a1, b1, a2, b2)
         double precision, intent(in) :: x1, a1, b1, a2, b2
      end function map
   end interface


   double precision :: x, y

   x = map(dble(rand()), 0d0, 1d0, a, b)
   y = cotasup * rand()

   do while (y > distr(x))
      x = map(dble(rand()), 0d0, 1d0, a, b)
      y = cotasup * rand()
   end do

   retval = x
end function generate_random_distributed_number

double precision function func_p(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   COMMON /CONSTS/PI
   double precision :: PI
   retval = 10d0/3d0 * dexp(-x) * dsin(x)**3 / (1+dexp(-PI))
end function func_p

function random_metropolis(dims, quantity, distr) result(retval)
   implicit none

   interface
      function rand_vect(dims1, a, b) result(retval1)
         integer, intent(in) :: dims1
         double precision, intent(in) :: a, b
         double precision, dimension(dims1) :: retval1
      end function rand_vect
      double precision function distr(xs)
         implicit none
         double precision, dimension(:), intent(in) :: xs
      end function distr
      function sum_vectors(dim, vec_a, vec_b)result(retval1)
         integer, intent(in) :: dim
         double precision, dimension(dim), intent(in) :: vec_a, vec_b
         double precision, dimension(dim) :: retval1
      end function sum_vectors
   end interface

   integer, intent(in) :: dims, quantity
   double precision, dimension(quantity,dims) :: retval
   double precision, dimension(dims) :: current_vector, new_vector
   integer :: i, k
   double precision :: r, p

   retval = 0
   current_vector = rand_vect(dims, -1d0, 1d0)

   do i = 2, quantity
      new_vector = rand_vect(dims, -1d0, 1d0)
      p = dble(rand())
      r = distr(new_vector)/distr(current_vector)
      if (r > p) then
         do k = 1, dims ! assign vector
            retval(i, k) = new_vector(k)
         end do
         current_vector = new_vector
      else
         do k = 1, dims ! assign vector
            retval(i, k) = retval(i-1, k)
         end do
      end if
   end do

end function random_metropolis

function rand_vect(dims, a, b) result(retval)
   implicit none

   integer, intent(in) :: dims
   double precision, intent(in) :: a, b
   double precision, dimension(dims) :: retval
   interface
      double precision function map(x1, a1, b1, a2, b2)
         double precision, intent(in) :: x1, b1, a1, a2, b2
      end function map
   end interface

   integer :: i

   do i = 1, dims
      retval(i) = map(dble(rand()), 0d0, 1d0, a, b)
   end do
end function rand_vect

double precision function integral_1(x) result(retval)
   implicit none

   COMMON /CONSTS/PI
   double precision :: PI

   double precision, intent(in) :: x

   retval = dsqrt(PI**2-x**2)
end function integral_1

double precision function integral_2(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   retval = ((x**2)*dsin(x) - x**3)*(dcos(x)**2)*dsin(x)
end function integral_2

double precision function integral_3(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   retval = dexp(-dabs(x)) * x**2*dsin(x)**2
end function integral_3

double precision function integral_4(x) result(retval)
   implicit none

   double precision, intent(in) :: x
   COMMON /CONSTS/PI
   double precision :: PI

   retval = dexp(-x**2/2) * dcos(x)**2*(PI + 4*x**2)
end function integral_4

double precision function integral_5(x) result(retval)
   implicit none

   double precision, intent(in) :: x
   COMMON /CONSTS/PI
   double precision :: PI

   retval = dexp(-x**2) * dsin(x)**2*x**2
end function integral_5

double precision function func_g(xs) result(retval)
   implicit none
   double precision, dimension(:), intent(in) :: xs
   retval = dexp(-(xs(1)**2+xs(2)**2+2*xs(3)**2+xs(4)**2+2*xs(5)**2))
end function func_g

double precision function integral_6(xs) result(retval)
   implicit none
   double precision, intent(in), dimension(:) :: xs
   COMMON /CONSTS/PI
   double precision :: PI
   retval = dexp(xs(2)*dcos(xs(2)+xs(3)-xs(5)))*(PI*xs(3)**2*xs(4)**2*xs(5)**3+dcos(xs(3)+xs(4)-2*xs(1))**2*xs(3)*dsin(xs(5)))
end function integral_6
