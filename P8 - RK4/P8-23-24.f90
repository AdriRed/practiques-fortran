program prepra8
   implicit none

   interface
      subroutine wave_function_calculator(energy, xini, xend, steps, yini1, normalize, xs1, phis1)
         double precision, intent(in) :: energy, xini, xend
         double precision, dimension(2), intent(in) :: yini1
         integer, intent(in) :: steps
         logical, intent(in) :: normalize

         double precision, dimension(steps), intent(out) :: xs1, phis1
      end subroutine wave_function_calculator
      subroutine eigenvalues_finder(E1, E2, eps, xini, xend, steps, yini1, write_eigenvalues, E3, xs1, phis1)
         implicit none

         integer, intent(in) :: steps
         double precision, intent(in) :: E1, E2, xini, xend, eps
         double precision, dimension(2), intent(in) :: yini1
         logical, intent(in) :: write_eigenvalues

         double precision, intent(out) :: E3
         double precision, intent(out), dimension(steps) :: xs1, phis1
      end subroutine eigenvalues_finder
      double precision function map(x, a, b, c, d) result(retval)
         implicit none
         double precision, intent(in) :: x, a, b, c, d

      end function map
   end interface

   common /potential/delta, V_0, beta, epsilon
   double precision :: delta, V_0, beta, epsilon

   common /constants/h2me
   double precision :: h2me

   double precision, dimension(2) :: yini
   double precision :: L, eigenvalue, probability, step
   double precision, dimension(400) :: xs, phis
   double precision, dimension(6) :: energies
   double precision, dimension(3) :: betas
   integer :: i, k, msigma_index, psigma_index


   delta = 0.4d0
   V_0 = -50d0
   L = 14d0
   h2me = 3.80995d0
   beta = 0d0
   epsilon = 1d0

   energies = (/-31d0, -30d0, -14d0, -13d0, -4d0, -3.5d0/)
   yini = (/0.0d0, 2d-6/)

   open(1, file="P8-23-24-res.dat")

   write(1, *) "# apartat 1"
   do i = 1, 4
      write(1, "(A, I1)") "#phi", i
      write(1, *) "# x, phi"
      call wave_function_calculator(energies(i), -L/2d0, L/2d0, 400, yini, .false., xs, phis)

      do k = 1, 400
         write(1, "(E20.12, E20.12)") xs(k), phis(k)
      end do

      write(1, *)
      write(1, *)
   end do

   write(1, *) "# apartat 2"

   do i = 1, 6, 2
      write(1, "(A, I1)") "#autovalors ", (i+1)/2
      call eigenvalues_finder(energies(i), energies(i+1), 1d-6, -L/2d0, L/2d0, 400, yini, .true., &
         eigenvalue, xs, phis)
      write(1, *)
      write(1, *)

      write(1, "(A, I1)") "#phi normalitzada ", (i+1)/2
      do k = 1, 400
         write(1, "(E20.12, E20.12)") xs(k), phis(k)
      end do
      write(1, *)
      write(1, *)

   end do

   betas = (/0d0, 5d0, 15d0/)

   write(1, *) "# apartat 3"

   msigma_index = int(map(-delta, -L/2d0, L/2d0, 1d0, 400d0))
   psigma_index = int(map(delta, -L/2d0, L/2d0, 1d0, 400d0))
   step = L/400d0

   open(2, file="P8-23-24-res1.dat")
   write(2, *) "# probabilitat"
   write(2, *) "# beta, prob"
   do i = 1, 3
      write(1, "(A, F5.1)") "# beta ", betas(i)
      write(1, *) "# x, phi"
      beta = betas(i)
      call eigenvalues_finder(energies(1), energies(2), 1d-6, -L/2d0, L/2d0, 400, yini, .false., &
         eigenvalue, xs, phis)

      probability = sum(phis(msigma_index:psigma_index)**2)*step

      write(2, "(F5.1, E20.12)") betas(i), probability

      do k = 1, 400
         write(1, "(E20.12, E20.12)") xs(k), phis(k)
      end do

      write(1, *)
      write(1, *)
   end do


   close(1)
   close(2)

end program prepra8

double precision function map(x, a, b, c, d) result(retval)
   implicit none
   double precision, intent(in) :: x, a, b, c, d

   retval = c + (x-a)*(d-c)/(b-a)

end function map

subroutine wave_function_calculator(energy, xini, xend, steps, yini, normalize, xs, phis)
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
            function edos(x, ys, nequs1) result(retval)
               implicit none

               integer, intent(in) :: nequs1
               double precision, intent(in) :: x
               double precision, dimension(nequs1), intent(in) :: ys
               double precision, dimension(nequs1) :: retval

            end function edos
         end interface
      end subroutine rungekutta4
      function schrodinger_equation(x, ys, nequs) result(retval)
         implicit none

         integer, intent(in) :: nequs
         double precision, intent(in) :: x
         double precision, dimension(nequs), intent(in) :: ys
         double precision, dimension(nequs) :: retval
      end function schrodinger_equation
   end interface

   double precision, intent(in) :: energy, xini, xend
   double precision, dimension(2), intent(in) :: yini
   integer, intent(in) :: steps
   logical, intent(in) :: normalize

   double precision, dimension(steps), intent(out) :: xs, phis

   double precision, dimension(2, steps) :: temp_result_ys

   common /schrodinger/E
   double precision :: E
   double precision :: normalization_constant, step

   step = (xend-xini)/dble(steps)

   E = energy
   call rungekutta4(steps, xini, xend, 2, yini, schrodinger_equation, xs, temp_result_ys)

   if (normalize) then
      ! sumatori de tots els valors de la funció pel pas
      normalization_constant = sum(temp_result_ys(1, 1:steps)**2)*step
      phis(1:steps) = temp_result_ys(1, 1:steps)/dsqrt(normalization_constant)
   else
      phis(1:steps) = temp_result_ys(1, 1:steps)
   end if

end subroutine wave_function_calculator

subroutine eigenvalues_finder(E1, E2, eps, xini, xend, steps, yini, write_eigenvalues, E3, xs, phis)
   implicit none

   integer, intent(in) :: steps
   double precision, intent(in) :: E1, E2, xini, xend, eps
   double precision, dimension(2), intent(in) :: yini
   logical, intent(in) :: write_eigenvalues

   double precision, intent(out) :: E3
   double precision, intent(out), dimension(steps) :: xs, phis

   interface
      subroutine rungekutta4(steps1, xini1, xend1, nequs, yini1, edos, result_xs, result_ys)
         implicit none

         integer, intent(in) :: nequs, steps1
         double precision, intent(in) :: xini1, xend1
         double precision, dimension(nequs), intent(in) :: yini1
         double precision, dimension(nequs,steps1), intent(out) :: result_ys
         double precision, dimension(steps1), intent(out) :: result_xs

         interface
            function edos(x, ys, nequs1) result(retval1)
               implicit none

               integer, intent(in) :: nequs1
               double precision, intent(in) :: x
               double precision, dimension(nequs1), intent(in) :: ys
               double precision, dimension(nequs1) :: retval1

            end function edos
         end interface
      end subroutine rungekutta4
      function schrodinger_equation(x, ys, nequs) result(retval1)
         implicit none

         integer, intent(in) :: nequs
         double precision, intent(in) :: x
         double precision, dimension(nequs), intent(in) :: ys
         double precision, dimension(nequs) :: retval1
      end function schrodinger_equation
   end interface

   COMMON /SCHRODINGER/E
   double precision :: E

   double precision :: phi1, phi2, phi3, current_E1, current_E2, normalization_constant, step
   double precision, dimension(2,steps) :: temp_result_ys
   double precision, dimension(steps) :: temp_result_xs
   integer :: i

   i = 1
   current_E1 = E1
   current_E2 = E2
   phi3 = eps+1d0
   step = (xend - xini)/steps

   do while (abs(phi3) > eps)
      E = current_E1
      call rungekutta4(steps, xini, xend, 2, yini, schrodinger_equation, temp_result_xs, temp_result_ys)
      phi1 = temp_result_ys(1, steps)

      E = current_E2
      call rungekutta4(steps, xini, xend, 2, yini, schrodinger_equation, temp_result_xs, temp_result_ys)
      phi2 = temp_result_ys(1, steps)

      E3 = (current_E1 * phi2 - current_E2 * phi1)/(phi2-phi1)
      E = E3
      call rungekutta4(steps, xini, xend, 2, yini, schrodinger_equation, temp_result_xs, temp_result_ys)
      phi3 = temp_result_ys(1, steps)

      if (write_eigenvalues) then
         write(1, "(I3, E20.12)") i, E3
      end if

      if (abs(phi3) > eps) then
         current_E1 = current_E2
         current_E2 = E3
         i = i + 1
      end if
   end do

   ! normalització
   normalization_constant = sum(temp_result_ys(1, 1:steps)**2)*step
   phis(1:steps) = temp_result_ys(1, 1:steps)/dsqrt(normalization_constant)
   xs(1:steps) = temp_result_xs(1:steps)


end subroutine eigenvalues_finder

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

!            { (1) -> dy/dx <- el valor de retval(1) es su derivada
!            | (2) -> d^2(y)/dx^2
!  retval =  | (3) -> d^3(y)/dx^3
!            | ...
!            { (nequs) -> d^(nequs)(y)/dx^(nequs)
function schrodinger_equation(x, ys, nequs) result(dy)
   implicit none

   integer, intent(in) :: nequs
   double precision, intent(in) :: x
   double precision, dimension(nequs), intent(in) :: ys
   double precision, dimension(nequs) :: dy

   interface
      double precision function potencial(x1)
         implicit none
         double precision, intent(in) :: x1
      end function potencial
   end interface

   common /SCHRODINGER/E
   double precision :: E

   common /constants/h2me
   double precision :: h2me

   dy(1) = ys(2)
   dy(2) = 1d0/h2me*(potencial(x)-E)*ys(1)

end function schrodinger_equation

double precision function potencial(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   common /potential/delta, V_0, beta, epsilon
   double precision :: delta, V_0, beta, epsilon
   retval = V_0*dsinh(2d0)/(dcosh(2d0)+dcosh(x/delta))+beta*dsin(x/epsilon)
end function potencial
