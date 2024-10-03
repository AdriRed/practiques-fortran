program practica10
   implicit none

   interface
      subroutine crank_nicholson(length11, position_steps1, time_limit1, time_steps1, &
         kcp, initial_conditions1, extra_contribution, unit, output1)
         implicit none

         interface
            double precision function extra_contribution(ix, t, length1, phi)
               implicit none
               integer, intent(in) :: length1, ix
               double precision, intent(in) :: t
               double precision, dimension(0:length1), intent(in) :: phi
            end function extra_contribution
         end interface
         double precision, intent(in) :: length11, time_limit1, kcp
         double precision, dimension(0:position_steps1), intent(in) :: initial_conditions1
         integer, intent(in) :: position_steps1, time_steps1
         integer, intent(in) :: unit
         double precision, dimension(0:position_steps1), intent(out) :: output1
      end subroutine crank_nicholson
      subroutine jacobi_method(length1, rho, h1, input_matrix, output_matrix)
         implicit none

         interface
            double precision function rho(i1, length11, input_matrix1)
               integer, intent(in) :: i1, length11
               double precision, intent(in), dimension(0:length11) :: input_matrix1
            end function
         end interface

         integer, intent(in) :: length1
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(0:length1) :: input_matrix
         double precision, intent(out), dimension(0:length1) :: output_matrix
      end subroutine jacobi_method
      subroutine poisson_equation_algorithm(method_used, length1, input_matrix, epsilon, &
         rho, h11, output_matrix)
         implicit none

         integer, intent(in) :: length1
         double precision, intent(in), dimension(0:length1) :: input_matrix
         double precision, intent(in) :: epsilon, h11
         double precision, intent(out), dimension(0:length1) :: output_matrix

         interface
            double precision function rho(i1, length11, input_matrix1)
               integer, intent(in) :: i1, length11
               double precision, intent(in), dimension(0:length11) :: input_matrix1
            end function
            subroutine method_used(length11, rho, h1, input_matrix1, output_matrix1)
               implicit none
               interface
                  double precision function rho(i1, length111, input_matrix11)
                     integer, intent(in) :: i1, length111
                     double precision, intent(in), dimension(0:length111) :: input_matrix11
                  end function
               end interface
               integer, intent(in) :: length11
               double precision, intent(in) :: h1
               double precision, intent(in), dimension(0:length11) :: input_matrix1
               double precision, intent(out), dimension(0:length11) :: output_matrix1
            end subroutine method_used
         end interface

      end subroutine poisson_equation_algorithm
      double precision function edo_right_side(i1, length1, input_matrix1)
         integer, intent(in) :: i1, length1
         double precision, intent(in), dimension(0:length1) :: input_matrix1
      end function edo_right_side
      subroutine write_crank_nicholson_step(unit, length1, time, position_step1, temperatures)
         implicit none

         integer, intent(in) :: length1, unit
         double precision, intent(in) :: time, position_step1
         double precision, dimension(0:length1), intent(in) :: temperatures
      end subroutine write_crank_nicholson_step
      double precision function extra_contribution(ix, t, length1, phi)
         implicit none

         integer, intent(in) :: length1, ix
         double precision, intent(in) :: t
         double precision, dimension(0:length1), intent(in) :: phi

      end function extra_contribution

   end interface

   common /cranknicholson/alpha, beta
   double precision :: alpha, beta

   double precision, dimension(3) :: betas
   double precision, dimension(0:50) :: initial_conditions, output

   double precision :: time_limit, length, ambient_temp, furnace_temp, x, alpha_fe, alpha_au
   integer :: time_steps, position_steps, i

   alpha_fe = 2.2d-5
   alpha_au = 1.29d-4

   length = 150
   time_limit = 4500d0
   furnace_temp = 280
   ambient_temp = 22
   time_steps = 6000
   position_steps = 50
   alpha = alpha_fe
   betas = (/0.00004, 0.0003, 0.0025/)


   open(1, file="P10-23-24-res-1.dat")
   do i = 1, 3
      beta = betas(i)
      write(1, *) "# temperature limit infinite, beta = ", beta
      initial_conditions = 0
      initial_conditions(position_steps) = furnace_temp-ambient_temp
      call poisson_equation_algorithm(jacobi_method, position_steps, initial_conditions, &
         1d-5, edo_right_side, length/position_steps, output)

      x = 0
      !                                               same as infinity
      call write_crank_nicholson_step(1, position_steps, -log(x), length/position_steps, output)
   end do
   close(1)

   open(2, file="P10-23-24-res-2a.dat")

   beta = 2d-4
   call crank_nicholson(length, position_steps, time_limit, time_steps, alpha, initial_conditions, &
      extra_contribution, 2, output)
   close(2)

   open(3, file="P12-23-24-res-2b.dat")
   alpha = alpha_fe
   beta = 15d-4

   write(3, *) "# alpha_fe"
   call crank_nicholson(length, position_steps, time_limit, time_steps, alpha, initial_conditions, &
      extra_contribution, 3, output)

   write(3, *) "# alpha_au"
   alpha = alpha_au
   call crank_nicholson(length, position_steps, time_limit, time_steps, alpha, initial_conditions, &
   extra_contribution, 3, output)
   close(3)

   open(4, file="P12-23-24-res-2c.dat")
   alpha = alpha_au
   beta = 2d-3

   call crank_nicholson(length, position_steps, time_limit, time_steps, alpha, initial_conditions, &
      extra_contribution, 4, output)
   close(3)

end program practica10

double precision function edo_right_side(i, length, input_matrix) result(retval)
   implicit none

   integer, intent(in) :: i, length
   double precision, intent(in), dimension(0:length) :: input_matrix
   common /cranknicholson/alpha, beta
   double precision :: alpha, beta


   retval = beta/alpha * input_matrix(i)
end function edo_right_side


subroutine poisson_equation_algorithm(method_used, length, input_matrix, epsilon, rho, h, output_matrix)
   implicit none

   interface
      double precision function rho(i1, length1, input_matrix1)
         integer, intent(in) :: i1, length1
         double precision, intent(in), dimension(0:length1) :: input_matrix1
      end function
      subroutine method_used(length1, rho, h1, input_matrix1, output_matrix1)
         implicit none
         interface
            double precision function rho(i1, length11, input_matrix11)
               integer, intent(in) :: i1, length11
               double precision, intent(in), dimension(0:length11) :: input_matrix11
            end function
         end interface
         integer, intent(in) :: length1
         double precision, intent(in) :: h1
         double precision, intent(in), dimension(0:length1) :: input_matrix1
         double precision, intent(out), dimension(0:length1) :: output_matrix1
      end subroutine method_used
   end interface

   integer, intent(in) :: length
   double precision, intent(in), dimension(0:length) :: input_matrix
   double precision, intent(in) :: epsilon, h
   double precision, intent(out), dimension(0:length) :: output_matrix

   double precision, dimension(0:length) :: current_matrix
   integer :: i

   call method_used(length, rho, h, input_matrix, output_matrix)
   output_matrix = 0
   i = 1

   current_matrix = input_matrix
   call method_used(length, rho, h, current_matrix, output_matrix)

   do while (maxval(dabs(current_matrix - output_matrix)) > epsilon)
      i = i+1
      current_matrix = output_matrix
      call method_used(length, rho, h, current_matrix, output_matrix)
   end do
end subroutine poisson_equation_algorithm

subroutine jacobi_method(length, rho, h, input_matrix, output_matrix)
   implicit none

   interface
      double precision function rho(i1, length1, input_matrix1)
         integer, intent(in) :: i1, length1
         double precision, intent(in), dimension(0:length1) :: input_matrix1
      end function
   end interface

   integer, intent(in) :: length
   double precision, intent(in) :: h
   double precision, intent(in), dimension(0:length) :: input_matrix
   double precision, intent(out), dimension(0:length) :: output_matrix

   integer :: i
   output_matrix = input_matrix

   output_matrix = input_matrix

   do i = 1, length-1
         output_matrix(i) = 0.5d0 * (input_matrix(i+1) + input_matrix(i-1) + &
            h**3*rho(i, length, input_matrix))
   end do
end subroutine jacobi_method


double precision function extra_contribution(ix, t, length1, phi) result(retval)
   implicit none

   integer, intent(in) :: length1, ix
   double precision, intent(in) :: t
   double precision, dimension(0:length1), intent(in) :: phi

   common /cranknicholson/alpha, beta
   double precision :: alpha, beta

   retval = -beta * phi(ix)
end function extra_contribution

subroutine crank_nicholson(length, position_steps, time_limit, time_steps, kcp, &
   initial_conditions, extra_contribution, unit, output)
   implicit none

   interface
      double precision function extra_contribution(ix, t, length1, phi)
         implicit none
         integer, intent(in) :: length1, ix
         double precision, intent(in) :: t
         double precision, dimension(0:length1), intent(in) :: phi
      end function extra_contribution
      subroutine crank_nicholson_step(current_temp, time_index, time_step1, r1, &
         a_maindiag1, a_lowerdiag1, a_upperdiag1, extra_contribution, length2, next_temp)
         implicit none

         interface
            double precision function extra_contribution(ix, t, length1, phi)
               implicit none
               integer, intent(in) :: length1, ix
               double precision, intent(in) :: t
               double precision, dimension(0:length1), intent(in) :: phi
            end function extra_contribution
         end interface

         integer, intent(in) :: length2, time_index
         double precision, dimension(0:length2), intent(in) :: current_temp
         double precision, dimension(1:length2-1), intent(in) :: a_maindiag1, a_lowerdiag1, a_upperdiag1
         double precision, intent(in) :: r1, time_step1

         double precision, dimension(0:length2), intent(out) :: next_temp
      end subroutine crank_nicholson_step
      subroutine write_crank_nicholson_step(unit, length1, time, position_step1, temperatures)
         implicit none

         integer, intent(in) :: unit, length1
         double precision, intent(in) :: time, position_step1
         double precision, dimension(0:length1), intent(in) :: temperatures
      end subroutine write_crank_nicholson_step
   end interface

   double precision, intent(in) :: length, time_limit, kcp
   integer, intent(in) :: position_steps, time_steps, unit
   double precision, dimension(0:position_steps), intent(in) :: initial_conditions
   double precision, dimension(0:position_steps), intent(out) :: output

   double precision, dimension(1:position_steps-1) :: a_maindiag, a_lowerdiag, a_upperdiag
   integer :: i

   double precision :: position_step, time_step, r
   double precision, dimension(0:position_steps) :: old_temperature

   position_step = length/position_steps
   time_step = time_limit/time_steps

   r = kcp*time_step/(position_step**2)

   ! definition of the A matrix
   a_upperdiag(1) = 0d0
   a_lowerdiag(position_steps-1) = 0d0

   a_maindiag = 2*(1d0+r)
   a_upperdiag(2:position_steps-1) = -r
   a_lowerdiag(1:position_steps-2) = -r

   old_temperature = initial_conditions

   if (unit > 0) then
      call write_crank_nicholson_step(unit, position_steps, 0d0, position_step, old_temperature)
      do i = 1, time_steps
         call crank_nicholson_step(old_temperature, i, time_step, r, a_maindiag, a_lowerdiag, &
            a_upperdiag, extra_contribution, position_steps, output)
         old_temperature = output
         call write_crank_nicholson_step(unit, position_steps, i*time_step, position_step, output)
      end do
   end if
end subroutine

subroutine write_crank_nicholson_step(unit, length, time, position_step, temperatures)
   implicit none

   integer, intent(in) :: length, unit
   double precision, intent(in) :: time, position_step
   double precision, dimension(0:length), intent(in) :: temperatures

   integer :: i

   write(unit, *) "# time = ", time
   write(unit, *) "# x, T-T_amb "
   do i = 0, length
      write(unit, "(E20.12, E20.12)") position_step*i, temperatures(i)
   end do
   write(unit, *)
   write(unit, *)


end subroutine write_crank_nicholson_step

! metodo sencillo de paso de crank-nicholson.
! extra contribution es S(x) en 6.46
! r = k*Δt/(c*ρ*(Δx)^2)
subroutine crank_nicholson_step(current_temp, time_index, time_step, r, a_maindiag, a_lowerdiag, &
   a_upperdiag, extra_contribution, position_steps, next_temp)
   implicit none

   interface
      subroutine tridiagonalization(length1, right_side_b, a_lower_diag, a_main_diag, a_upper_diag, phi)
         integer, intent(in) :: length1
         double precision, dimension(0:length1-1), intent(in) :: right_side_b, a_lower_diag, &
            a_main_diag, a_upper_diag
         double precision, dimension(0:length1), intent(inout) :: phi
      end subroutine tridiagonalization

      double precision function extra_contribution(ix, t, length1, phi)
         implicit none
         integer, intent(in) :: length1, ix
         double precision, intent(in) :: t
         double precision, dimension(0:length1), intent(in) :: phi
      end function extra_contribution
   end interface

   integer, intent(in) :: position_steps, time_index
   double precision, dimension(0:position_steps), intent(in) :: current_temp
   double precision, dimension(1:position_steps-1), intent(in) :: a_maindiag, a_lowerdiag, a_upperdiag
   double precision, intent(in) :: r, time_step

   double precision, dimension(0:position_steps), intent(out) :: next_temp
   double precision, dimension(1:position_steps-1) :: b_diag
   double precision :: time
   integer :: i

   time = time_index*time_step
   ! calculamos la contribucion espacial de los valores actuales en b_diag (valor derecho de la igualdad 6.56 de los apuntes)
   do i = 1, position_steps-1
      b_diag(i) = ( 2*current_temp(i)*(1d0-r) + r*current_temp(i+1) + r*current_temp(i-1)) &
         + extra_contribution(i, time, position_steps, current_temp) * time_step
   end do

   ! condiciones de contorno
   b_diag(1) = b_diag(1) + 2*r*current_temp(0)
   b_diag(position_steps-1) = b_diag(position_steps-1) + 2*r*current_temp(position_steps)

   next_temp = current_temp

   call tridiagonalization(position_steps, b_diag, a_lowerdiag, a_maindiag, a_upperdiag, next_temp)

end subroutine crank_nicholson_step

subroutine tridiagonalization(length, right_side_b, a_lower_diag, a_main_diag, a_upper_diag, phi)
   implicit none

   integer, intent(in) :: length
   double precision, dimension(1:length-1), intent(in) :: right_side_b, a_lower_diag, a_main_diag, a_upper_diag
   double precision, dimension(0:length), intent(inout) :: phi ! intent(inout) para obtener las condiciones de contorno (extremo del vector)

   double precision, dimension(1:length-1) :: alphas, betas, gammas
   integer :: i

   ! asignamos los valores iniciales y precalculados: texto debajo de 6.64 y 6.63
   alphas(1) = 0
   betas(1) = phi(0)

   ! expresiones 6.63 y 6.64
   do i = 1, length-2
      gammas(i) = -1d0/(a_main_diag(i) + a_lower_diag(i)*alphas(i))
      
      alphas(i+1) = gammas(i)*a_upper_diag(i)
      betas(i+1) = gammas(i)*(a_lower_diag(i)*betas(i)-right_side_b(i))
   end do

   ! expresion 6.62 con expresiones 6.64
   do i = 1, length-2
      phi(i+1) = alphas(i) *phi(i) + betas(i)
   end do

end subroutine
