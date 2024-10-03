! PROGRAMA SUMA
! BJD SEP 2015
!
! --- A,B NÚMEROS INTRODUCIDOS POR LA CONSOLA
!
    IMPLICIT NONE
PROGRAM main
    integer :: result

    result = mult(3, 6)
    
    print * result

END PROGRAM main

integer function mult(a, b) result(retval)
    integer :: a, b
    retval = a*b
    return
end function mult