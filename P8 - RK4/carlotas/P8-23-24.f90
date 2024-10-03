!CARLOTA LEOPOLD
!DIMARTS
PROGRAM passat
    IMPLICIT NONE
    DOUBLE PRECISION::V0, delta, alpha, C, phi0, dphi0, L, xi, xf, x, dx, res_phi, integral, nequs
    DOUBLE PRECISION:: E, E1, E2, E3, E4, E5, E6, Ef, Ef1, Ef2, Ef3, error
    INTEGER::N, i
    DOUBLE PRECISION, ALLOCATABLE::phi_1(:), phi_2(:), phi_3(:), phi_4(:)
    DOUBLE PRECISION, EXTERNAL:: V
    COMMON/POTENCIAL/V0, delta, alpha
    COMMON/ENERGIA/E, C, phi0, dphi0

    OPEN(222, file="P8-23-24-res.dat")

    C=7.6199d0 !eV*A**2
    delta=0.05d0 !A
    V0=-20.d0 !eV
    
    phi0=0.d0 !A**-1/2
    dphi0=2.d-6 !A**-3/2

    L=8.d0 !A
    xi=-L/2
    xf=L/2

    N=400
    nequs=2

    !Apartat 1
    ALLOCATE(phi_1(N), phi_2(N), phi_3(N), phi_4(N))
    E1=-21.d0 !eV
    E2=-20.5d0 !eV
    E3=-14.d0 !eV
    E4=-13.d0 !eV

    CALL Integrar(xi, xf, N, E1, phi_1, V, res_phi, integral)
    CALL Integrar(xi, xf, N, E2, phi_2, V, res_phi, integral)
    CALL Integrar(xi, xf, N, E3, phi_3, V, res_phi, integral)
    CALL Integrar(xi, xf, N, E4, phi_4, V, res_phi, integral)

    dx = (xf-xi)/dble(N)
    WRITE(222, *) "#x,      phi_1(i),       phi_2(i),       phi_3(i),       phi_4(i)"
    DO i=1, N 
        x = xi + dx*i
        WRITE(222, *) x, phi_1(i), phi_2(i), phi_3(i), phi_4(i)
    END DO

    WRITE(222, *) ""
    WRITE(222, *) ""

    !APARTAT 2
    E5=-8.d0 !eV
    E6=-7.5d0 !eV
    error = 10.d-6 !A^(-1/2)

    CALL  Tir(xi, xf, E1, E2, Ef1, N, error, V)
    CALL  Tir(xi, xf, E3, E4, Ef2, N, error, V)
    CALL  Tir(xi, xf, E5, E6, Ef3, N, error, V)
    
    CLOSE(222) 
END PROGRAM

!SUBRUTINES--------------------------------------------------------------------------
!Subrutina Runge Kutta 4
SUBROUTINE RungeKutta4order(x0, dx, funcin, dfuncout, nequs, funcio)
    IMPLICIT NONE
    DOUBLE PRECISION:: x0, dx, x2, x3, x4, E, C, phi0, dphi0
    INTEGER::NEQUS, i
    DOUBLE PRECISION, DIMENSION(NEQUS)::funcin, dfuncout, k1, y, k2, k3, k4
    DOUBLE PRECISION, EXTERNAL::funcio
    COMMON/ENERGIA/E, C, phi0, dphi0

    !k1=f(x0, y0)
    CALL EDO(NEQUS, x0, funcin, k1, funcio)

    !k2
    y = funcin + (dx/2.d0)*k1
    x2 = x0 + dx/2.d0
    CALL EDO(NEQUS, x2, y, k2, funcio)

    !k3
    y = funcin + (dx/2.d0)*k2
    x3 = x0 + dx/2.d0
    CALL EDO(NEQUS, x3, y, k3, funcio)

    !k4
    y = funcin + dx*k3
    x4 = x0 + dx
    CALL EDO(NEQUS, x4, y, k4, funcio)

    dfuncout = funcin + (dx/6.d0)*(k1 + 2.d0*k2 + 2.d0*k3 + k4)

RETURN
END

!Subrutina EDO
SUBROUTINE EDO(nequ, x, yinput, dyoutput, funcio)
    IMPLICIT NONE
    DOUBLE PRECISION::x, x0, yinput(nequ), dyoutput(nequ), E, C, phi0, dphi0
    INTEGER::nequ
    DOUBLE PRECISION, EXTERNAL::funcio
    COMMON/ENERGIA/E, C, phi0, dphi0

    !Pel cas Eq d'Schrodinger on y=(phi, dphi/dx)
    dyoutput(1) = yinput(2)
    dyoutput(2) = 2.d0*(funcio(x0) - E)*yinput(1)/C
RETURN
END

!Subrutina MÈTODE TRAPEZIS
!fent servir la regla de Trapezis amb 3k intervals.
SUBROUTINE trapezoidalrule(x1, x2, k, func, resultat)
    IMPLICIT NONE
    DOUBLE PRECISION::x1, x2, h, xi, resultat, func(k)
    INTEGER:: k, i, interval

    interval = 3**k
    h = (x2-x1)/interval
    resultat = (h/2.d0)*(func(1) + func(k))

    DO i = 1, interval-1
        xi = x1 + i*h
        resultat = resultat + h*(func(i+1))
    END DO
RETURN
END SUBROUTINE

!Subrutina Integrar
SUBROUTINE Integrar(xi, xf, N, E0, list_phi, funcio, res_phi, integral)
    IMPLICIT NONE
    DOUBLE PRECISION ::a, b, xi, xf, E, C, phi0, dphi0,  dx, integral, E0, res_phi, x, resultat
    DOUBLE PRECISION, DIMENSION(2) ::phi, phi_out
    DOUBLE PRECISION, DIMENSION(N) ::list_phi, list_phi2
    INTEGER ::N, i, j, nequs
    DOUBLE PRECISION, EXTERNAL::funcio
    COMMON/ENERGIA/E, C, phi0, dphi0
    a=xi
    b=xf

    !Valors inicials de la funcio i la derivada
    phi(1) = phi0
    phi(2) = dphi0 
    E=E0
    dx = (b-a)/dble(N) 
    NEQUS=2

    DO  i = 1, N !següents valors de phi 1
        x = a + dx*i !intervals per la variable x
        list_phi(i) = phi(1) 
        !list_phi2(i) = phi(1)**2
        CALL RungeKutta4order(x, dx, phi, phi_out, nequs, funcio)
        DO j = 1, NEQUS
            phi(j) = phi_out(j) !valors finals de la funcio i la derivada despres de rk4
        END DO
    END DO

    res_phi = phi_out(1) 
    
    !Integrem metode extern phil
    !CALL simpsontresvuit(xi, xf, N, list_phi, resultat)
    CALL trapezoidalrule(xi, xf, N, list_phi, resultat)
    resultat=integral

    !REVISAR PQ NOSE SI SHA DE FER SERVIR EL QUADRAT O NO

END SUBROUTINE

!NO ES MEU, NO EL FAIG SERVIR
!subrutina externa per calcular l'integral,ara enlloc de tenir una funcio externa a integrar
!tindrem una llista de valors
SUBROUTINE simpsontresvuit(xi, xf, k, funcio, resultat)
    IMPLICIT NONE
    !DEFINIM LES VARIABLES
    INTEGER:: m, k, interval
    DOUBLE PRECISION:: h, xi, xf, resultat, funcio(k)
    
    !NUM MAX INTERVALS
    interval = k
    
    !DISTANCIA ENTRE PUNTS
    h = ABS(xf-xi)/dble(interval)
    !definim els els extrems i els sumem a l'integral
    resultat = (funcio(1) + funcio(k))*h*(3.d0/8.d0)

    !els altres punts
    DO m = 0, interval-1 
        !PELS PUNTS QUE SIGUIN DIVISIBLES PER 3 MULTIPLICAREM 2X3/8
        IF (mod(m,3) .EQ. 0) THEN
                resultat = resultat + funcio(m+1)*h*(6.d0/8.d0)
        ELSE ! I ELS QUE NO SON NI DIVISIBLES PER 3 NI EXTREMS
                resultat = resultat + funcio(m+1)*h*(9.d0/8.d0)
        END IF
    END DO
RETURN
END SUBROUTINE 

!Subrutina Mètode Tir
SUBROUTINE Tir(xi, xf, E1, E2, Ef, N, error, funcio)
    IMPLICIT NONE
    DOUBLE PRECISION :: E1, E2, Ef, error, xi, xf, res_phi1, res_phi2, res_phi3
    DOUBLE PRECISION :: integral1, integral2, integral3, list_phi(N)
    INTEGER :: N, k, i
    DOUBLE PRECISION, EXTERNAL:: funcio

    xi = 0.d0
    xf = 1.d0
 
    DO i = 1, k
        !Integrem 
        CALL integrar(xi, xf, k, E1, list_phi, funcio, res_phi1, integral1)
        CALL integrar(xi, xf, k, E2, list_phi, funcio, res_phi2, integral2)

        ! definim E3 amb els resultats de la subrutina  
        Ef = (E1*res_phi2 - E2*res_phi1)/(res_phi2 - res_phi1)

        !integrem el valor d'Ef
        CALL integrar(xi, xf, k, Ef, list_phi, funcio, res_phi3, integral3)
        !escribim per cada iteracio el valor de E3 en el file (des d'aqui perque és més fàcil), fins que convergeixi
        WRITE(222,*) i, Ef

        !volem que convergeixi si el valor del resultat de rungekutta es menor que l'error
        IF (ABS(res_phi3) .LT. error) THEN
            WRITE(222,*) "Ha convergit", i, "Energia (eV) :", Ef
            EXIT
        ELSE
            !si no convergeix redefinim
            E1 = E2
            E2 = Ef
        END IF
    END DO
   
END SUBROUTINE

!FUNCIONS----------------------------------------------------------------------------
!Funció Potencial(x)
DOUBLE PRECISION FUNCTION V(x)
    IMPLICIT NONE
    DOUBLE PRECISION::x
    COMMON/POTENCIAL/V0, delta, alpha
    double precision :: V0, delta, alpha
    V= V0 * (DSINH(alpha/delta) /(DCOSH(alpha/delta)+DCOSH(x/delta)))
END FUNCTION