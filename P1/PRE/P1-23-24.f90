program PREPRA1
   implicit none

   common /ENERGY/E_1
   real :: E_1 = 3.72, calc_energy, calc_fermi_energy !eV
   integer :: N = 40, read_k_val
   integer :: k, i
   real :: e_k, e_fermi

   k = read_k_val()

   e_k = calc_energy(k)
   e_fermi = calc_fermi_energy(N)

   write (*, "(a, I2, a, F0.2, a)") "L'energia d'una particula amb una k = ",  k, " es: ", e_k, " eV"
   write (*, "(a, F0.2, a)") "L'energia de Fermi per N = 40 es: ", e_fermi, " eV"

   open(1, file='P1-23-24-res1.dat')

   do i = 1, 40
      e_fermi = calc_fermi_energy(i)
      write(1, *) i, e_fermi, (calc_fermi_energy(i*2)/e_fermi)
   end do
   read(*, *)

end program PREPRA1

real function calc_energy(k) result(energy)
   implicit none

   common /ENERGY/E_1
   real :: E_1
   integer, intent(in) :: k

   energy = E_1*k**2
end function calc_energy

real function calc_fermi_energy(n) result(fermi_energy)
   implicit none

   real :: calc_energy
   integer, intent(in) :: n
   integer k

   fermi_energy = 0
   do k = 1, N
      fermi_energy = fermi_energy + calc_energy(k)
   end do

end function calc_fermi_energy

integer function read_k_val() result(k)
   implicit none

   k = 0

   do while (k < 2 .OR. 40 < k)
      write (*, *) "Escriu un valor de k (min 2, max 40)"
      read (*, *) k
      if (k < 2 .OR. 40 < k) then
         write (*, "(a, I3.0, a)") "Valor", k, " es invalid"
      end if
   end do

end function read_k_val
