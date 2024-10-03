program test
   implicit none
   real :: test_bisect, bisect

   write(*, *) bisect(test_bisect, 0, 2.5, 0.0001)
   


end program test

real function test_bisect(x) result(retval)
   real, intent(in) :: x
   retval = sin(x)**2*x-1
end function test_bisect

real function bisect(f, a, b, epsilon) result(retval)
   real function :: f, a, b, epsilon
   real :: c, f_c
   integer :: maxiter

   maxiter = nint(log((b-a)/epsilon/log(2.))) + 1
   do index = 1, maxiter
      c = (a+b)/2.
      f_c = f(c)
      if (f_c == 0.) then
         error stop
      end if

      if (f_c * f(a) < 0.) then
         b = c
      else
         a = c
      endif

      if (b-a < epsilon) then
         c = b-a
         stop
      end if
   end do
   retval = c
end function bisect
