program prepra8
   implicit none

   interface
      double precision function eigenvalues_finder(E1, E2, eps, xini, xend, steps, yini) result(retval)
         implicit none

         integer, intent(in) :: steps
         double precision, intent(in) :: E1, E2, xini, xend, eps
         double precision, dimension(2), intent(in) :: yini
         double precision, dimension(2,steps) :: temp_results
      end function eigenvalues_finder
   end interface

   double precision, dimension(2) :: yini
   double precision :: En
   integer :: i

   yini = (/0.0d0, 0.15d0/)

   open(1, file="data.dat")

   do i = 1, 4
      write(1, "(A,1X,I2,1X,A)") "# nivell ", i, ", 20 pasos"
      En = eigenvalues_finder(4.9d0*i**2-2.4, 5d0*i**2-2.4, 1d-5, 0d0, 1d0, 20, yini)
      write(1, *)
      write(1, *)
      write(1, *)
      write(1, *)
      write(1, "(A,1X,I2,1X,A)") "# nivell ", i, " 400 pasos"
      En = eigenvalues_finder(4.9d0*i**2-2.4, 5d0*i**2-2.4, 1d-5, 0d0, 1d0, 400, yini)
      write(1, *)
      write(1, *)
      write(1, *)
      write(1, *)
   end do 

   close(1)

end program prepra8

double precision function eigenvalues_finder(E1, E2, eps, xini, xend, steps, yini) result(retval)
   implicit none

   integer, intent(in) :: steps
   double precision, intent(in) :: E1, E2, xini, xend, eps
   double precision, dimension(2), intent(in) :: yini

   interface
      subroutine rungekutta4(steps, xini, xend, nequs, yini, edos, result_xs, result_ys)
         implicit none

         integer, intent(in) :: nequs, steps
         double precision, intent(in) :: xini, xend
         double precision, dimension(nequs), intent(in) :: yini
         double precision, dimension(nequs,steps), intent(out) :: result_ys
         double precision, dimension(steps), intent(out) :: result_xs

         interface
            function edos(x, ys, nequs) result(retval)
               implicit none

               integer, intent(in) :: nequs
               double precision, intent(in) :: x
               double precision, dimension(nequs), intent(in) :: ys
               double precision, dimension(nequs) :: retval

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
      double precision function simpsonstresvuit(x1, x2, steps, ys) result(retval)
         implicit none

         double precision, intent(in) :: x1, x2
         double precision, dimension(steps), intent(in) :: ys
         integer, intent(in) :: steps
      end function simpsonstresvuit
      double precision function potencial(x) 
         double precision, intent(in) :: x 
      end function potencial
   end interface

   COMMON /SCHRODINGER/E
   double precision :: E

   double precision :: phi1, phi2, phi3, E3, current_E1, current_E2, sum_value, step
   double precision, dimension(2,steps) :: temp_result_ys
   double precision, dimension(steps) :: temp_result_xs
   integer :: i

   i = 1
   current_E1 = E1
   current_E2 = E2
   phi3 = eps+1d0
   step = (xend - xini)/steps

   write(1, *) "# Autovalors"
   write(1, *) "# i, E"
   
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
      
      ! codigo sucio     
      write(*, *) i, E3, phi3
      write(1, "(I3, E20.12)") i, E3

      if (abs(phi3) > eps) then
         current_E1 = current_E2
         current_E2 = E3
         i = i + 1
      end if
   end do

   sum_value = sum(temp_result_ys(1, 1:steps)**2)*step
   temp_result_ys(1, 1:steps) = temp_result_ys(1, 1:steps)/dsqrt(sum_value)
   write(*, *) "Phi found at ", phi3, " with E = ", E3

   write(1, *)
   write(1, *)

   write(1, *) "# autovectors"
   write(1, *) "# x, phi, phi^2"
   
   do i = 1, steps
      write(1, "(E20.12, E20.12, E20.12)") temp_result_xs(i), temp_result_ys(1, i), temp_result_ys(1, i)
   end do

   write(1, *)
   write(1, *)

   retval = E3
end function eigenvalues_finder

subroutine rungekutta4(steps, xini, xend, nequs, yini, edos, result_xs, result_ys)
   implicit none

   integer, intent(in) :: nequs, steps
   double precision, intent(in) :: xini, xend
   double precision, dimension(nequs), intent(in) :: yini
   double precision, dimension(nequs,steps), intent(out) :: result_ys
   double precision, dimension(steps), intent(out) :: result_xs

   interface
      function edos(x, ys, nequs) result(retval)
         implicit none

         integer, intent(in) :: nequs
         double precision, intent(in) :: x
         double precision, dimension(nequs), intent(in) :: ys
         double precision, dimension(nequs) :: retval

      end function edos
      subroutine rungekutta4step(x, h, funcin, dfuncout, nequs, edofuncio)
         implicit none
         integer, intent(in) :: nequs
         double precision, dimension(nequs), intent(in) :: funcin
         double precision, dimension(nequs), intent(out) :: dfuncout
         double precision, intent(in) :: x, h

         interface
            function edofuncio(x, ys, nequs) result(retval)
               implicit none

               integer, intent(in) :: nequs
               double precision, intent(in) :: x
               double precision, dimension(nequs), intent(in) :: ys
               double precision, dimension(nequs) :: retval

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
      function edofuncio(x, ys, nequs) result(retval)
         implicit none

         integer, intent(in) :: nequs
         double precision, intent(in) :: x
         double precision, dimension(nequs), intent(in) :: ys
         double precision, dimension(nequs) :: retval

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
      double precision function potencial(x)
         implicit none
         double precision, intent(in) :: x
      end function potencial
   end interface

   COMMON /SCHRODINGER/E
   double precision :: E

   dy(1) = ys(2)
   dy(2) = 2d0*(potencial(x)-E)*ys(1)

end function schrodinger_equation

double precision function potencial(x) result(retval)
   implicit none
   double precision, intent(in) :: x
   retval = -2.4d0
end function potencial
