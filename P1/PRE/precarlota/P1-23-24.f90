! PREPRA 1
! Carlota Leopold

program prepra1
    implicit none

    !variables
    INTEGER :: k, N
    REAL :: E_k
    REAL :: E_F
    !funcions
    REAL :: energia_fermi

!Apartat 1
    READ(*, *) k

    IF (2<=k .AND. k<= 40) THEN
        E_k=3.72*k**2
        WRITE (*, *) E_k
    ELSE 
        WRITE (*, *) "Aquest número no està entre 2 i 40. "

    END IF

!Apartat 2
    E_F = energia_fermi(40)
    WRITE(*, *) E_F

!Apartat 3
    OPEN(1, file="P1-23-24-res1.dat")

    do N = 1, 40
        WRITE(1, *) N, energia_fermi(N)
    end do

!Apartat 4
    OPEN(2, file="P1-23-24-res2.dat")

    do N = 1, 20
        WRITE(2, *) N, (energia_fermi(2*N)/energia_fermi(N)), 8.-6./N+6./N**2.
    end do

    read(*, *)
end program prepra1

real function energia_fermi(N) result(ret)
    implicit none
    integer, intent(in) :: N
    integer :: k
    REAL :: E_Fermi = 0
    E_Fermi = 0
    do k=1, N
        E_Fermi = E_Fermi + 3.72*k**2
    end do

    ret = E_Fermi
end function energia_fermi

