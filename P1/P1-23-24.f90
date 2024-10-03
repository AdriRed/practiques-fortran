program PRACTICA1
   implicit none

   COMMON /CONSTS/E

   real(kind = 8) :: E = 2.718281828, precalc, p_k, calc_s_m, calc_s_asim
   integer :: k, read_k_val, n

   k = read_k_val() 
   write(*, *) p_k(k) ! apartado 1
   write(*, *) calc_s_m(28, 65) !apartado 2

   open(1, file='P1-23-24-res1.dat')
   do n = 11, 331, 3
      precalc = calc_s_m(8, n)
      write(1, *) n, precalc, calc_s_asim(n), precalc / calc_s_asim(n) !apartados 3 y 4
   end do

   read (*, *) ! para mantener terminal abierto pq se me cierra al acabar xd
end program PRACTICA1

real(kind = 8) function calc_s_m(m, n) result(retval)
   implicit none

   real(kind = 8) :: p_k
   integer, intent(in) :: n, m
   integer k

   retval = 0
   do k = m, n
      retval = retval + p_k(k)
   end do

end function calc_s_m

real(kind = 8) function calc_s_asim(n) result(retval)
   implicit none
   integer, intent(in) :: n
   retval = 1./5. * n ** 3
end function calc_s_asim

real(kind = 8) function p_k(k) result(retval)
   integer, intent(in) :: k
   COMMON /CONSTS/E

   real(kind = 8) :: E
   
   retval = 3./5.*k**2 + E + 10*k
end function p_k

integer function read_k_val() result(k)
   implicit none

   k = 0

   do while (k < 15 .OR. 221 < k)
      write (*, *) "Escriu un valor de k (min 15, max 221)"
      read (*, *) k
      if (k < 15 .OR. 221 < k) then
         write (*, "(a, I3.0, a)") "Valor", k, " es invalid"
      end if
   end do

end function read_k_val
