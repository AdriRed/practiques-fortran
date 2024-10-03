program practca6
   implicit none

   interface
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

      double precision function quark_up(x)
         implicit none
         double precision, intent(in) :: x
      end function quark_up
      double precision function quark_down(x)
         implicit none
         double precision, intent(in) :: x
      end function quark_down
      double precision function atom_probability(x)
         implicit none
         double precision, intent(in) :: x
      end function atom_probability
      double precision function func_g(x)
         implicit none
         double precision, intent(in) :: x
      end function func_g
      double precision function mitja_circumferencia(x)
         implicit none
         double precision, intent(in) :: x
      end function mitja_circumferencia

      subroutine encert_error_integral(func, count, a, b, cotasup, integral,  error)
         implicit none

         integer, intent(in) :: count
         double precision, intent(in) :: a, b, cotasup
         double precision, intent(out) :: integral, error

         interface
            double precision function func(x1)
               double precision, intent(in) :: x1
            end function func
            double precision function map(x1, a1, b1, a2, b2) result(retval)
               double precision, intent(in) :: x1, b1, a1, a2, b2
            end function map
         end interface
      end subroutine encert_error_integral

   end interface

   COMMON /CONSTS/PI
   double precision :: PI = datan(1d0)*4d0
   COMMON /ATOM/L
   double precision :: L

   integer :: i
   double precision :: integral_quark_up, integral_quark_down
   double precision :: error_quark_up, error_quark_down

   double precision :: integral, error
   double precision, dimension(1000000) :: atom_positions

   L = PI

   call srand(20401021)

   open(1, file="P6-23-24-res.dat")

   write(1, *) "# montecarlo cru quarks up & down"
   write(1, *) "# N, I_u, err_u, I_d, err_d"

   write(*, *) "Montecarlo Cru"
   do i = 150, 45000, 150
      call montecarlocru(quark_up, i, 0d0, 1d0, integral_quark_up, error_quark_up)
      call montecarlocru(quark_down, i, 0d0, 1d0, integral_quark_down, error_quark_down)
      write(*, *) i
      write(1, "(I6, E20.12, E20.12, E20.12, E20.12, E20.12, E20.12)") i, &
         integral_quark_up, error_quark_up, &
         integral_quark_down, error_quark_down
   end do

   write(*, *) "Montecarlo Cru OK"
   write(1, *)
   write(1, *)

   do i = 1, 1000000
      atom_positions(i) = generate_random_distributed_number(atom_probability, -L, L, 1/L)
   end do

   write(1, *) "# montecarlo sample I_3"
   write(1, *) "# N, I_2, err_2"

   write(*, *) "Montecarlo Sample"
   call montecarlosample(func_g, atom_probability, atom_positions, 1000000, integral, error)

   write(*, *) "Montecarlo Sample OK"
   write(1, *)
   write(1, *)

   write(1, *) "# numero pi"
   write(1, *) "# N, I_3, err_3"

   write(*, *) "Pi"
   do i = 10000, 300000, 10000
      call encert_error_integral(mitja_circumferencia, i, -1d0, 1d0, 1d0, integral, error)
      integral = integral * 2
      write(*, *) i
      write(1, "(I6, E20.12, E20.12)") i, &
         integral, error
   end do
   write(*, *) "Pi OK"
   write(1, *)
   write(1, *)




   close(1)

end program practca6

! funcions de la pracitca

double precision function mitja_circumferencia(x) result(retval)
   implicit none
   double precision, intent(in) :: x

   retval = dsqrt(1-x**2)

end function mitja_circumferencia

double precision function quark_up(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   retval = 5.109d0 * x**(0.8002d0-1d0)*(1-x)**3
end function quark_up

double precision function quark_down(x) result(retval)
   implicit none

   double precision, intent(in) :: x

   retval = 3.058d0 * x**(0.803d0-1d0)*(1-x)**4
end function quark_down

double precision function atom_probability(x) result(retval)
   implicit none

   COMMON /CONSTS/PI
   double precision :: PI
   common /ATOM/L
   double precision :: L

   double precision, intent(in) :: x
   retval = 1/L*dsin(PI*(x-L)/(2d0*L))**2

end function atom_probability

double precision function func_g(x) result(retval)
   double precision, intent(in) :: x
   COMMON /CONSTS/PI
   double precision :: PI
   common /ATOM/L
   double precision :: L

   retval = dsin(8d0*PI*(x-L)/(2d0*L))**2

end function func_g

subroutine encert_error_integral(funcio, total, a, b, cotasup, integral, error)
   implicit none

   integer, intent(in) :: total
   double precision, intent(in) :: a, b, cotasup
   double precision, intent(out) ::  integral, error

   interface
      double precision function funcio(x)
         implicit none
         double precision, intent(in) :: x
      end function funcio
      double precision function map(x, a1, b1, a2, b2) result(retval)
         implicit none
         double precision, intent(in) :: x, b1, a1, a2, b2
      end function map
   end interface

   integer :: i, count
   double precision :: x, y
   double precision :: dcount, dtotal

   count = 0

   do i = 1, total
      x = map(dble(rand()), 0d0, 1d0, a, b)
      y = rand() * cotasup
      
      if (funcio(x) >= y) then
         count = count + 1
      end if
   end do

   dtotal = dble(total)
   dcount = dble(count)
   integral = cotasup * (b-a) * dcount/dtotal
   error = cotasup * (b-a)/dsqrt(dtotal)* dsqrt(dcount/dtotal*(1-dcount/dtotal))
end subroutine encert_error_integral

! funcions de la prepra


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
      if (modulo(i, 10000) == 0) then
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

double precision function map(x, a1, b1, a2, b2) result(retval)
   implicit none
   double precision, intent(in) :: x, b1, a1, a2, b2
   retval = (x-a1)/(b1-a1) * (b2-a2) + a2
end function map
