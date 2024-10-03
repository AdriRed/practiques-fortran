!PREPRA 3
!CARLOTA LEOPOLD

program PREPRA
   IMPLICIT none

   DOUBLE PRECISION:: xini, A, B, FA, FB, FC, valordf, preci, valorarrel, DIFF
   INTEGER:: nitera, maxiter
   EXTERNAL:: fun, biseccio, newtonraphson

   CALL biseccio(fun,  A, B, preci, nitera, valorarrel)
   CALL newtonraphson(fun, xini, preci, nitera, valorarrel)
END PROGRAM

!SUBRUTINA FUN
SUBROUTINE fun(x, valorf, valordf)
   IMPLICIT none
   DOUBLE PRECISION:: x, valorf, valordf, pi, P

   pi=4*ATAN(1.d0)
   P=-(57.*pi)/160.+(57./80.+(17.*pi)/20.)*x-(17./10.+pi/2.)*x**2+x**3

   valorf = SINH(x)*P
   valordf=COSH(x)*P+(3*x**2+(-5*pi-17.)*x/5+(57.+pi*68.)/80.)*SINH(x)

   RETURN
END

!SUBRUTINA NEWTONRAPHSON
SUBROUTINE newtonraphson(fun, xini, preci, nitera, valorarrel)
   IMPLICIT none
   DOUBLE PRECISION:: preci, delta
   DOUBLE PRECISION::xini, fxini, fpxini, valorarrel
   INTEGER:: nitera, maxiter
   EXTERNAL:: fun

   preci=0.0001
   maxiter=100

   xini=0.2d0

   OPEN(34, file='P3-23-24-res1.dat')
   print*,'hola'
   DO nitera=1, maxiter

      CALL fun(xini, fxini, fpxini)

      valorarrel=xini-fxini/fpxini
      delta=ABS(fxini/fpxini)

      WRITE(34, "(A35,I3,A35,f20.12)") 'Iteració', nitera, 'Valor x=', valorarrel

      IF (delta .le. preci) THEN
         WRITE(34, "(A35,f20.12)") "Precisió conseguida número d'iteracions=", delta
         WRITE(34, "(A35,f20.12)") "L'arrel és x=", valorarrel
         STOP
      END IF
      xini=valorarrel
   END DO

   WRITE(34, "(A35)") 'El problema no ha convergit.'
   CLOSE(34)
   RETURN
END

!SUBRUTINA BISECCIO
SUBROUTINE biseccio(fun, A, B, preci, nitera, valorarrel)
   IMPLICIT none
   DOUBLE PRECISION:: A, B, FA, FB, FC, valordf, preci, valorarrel, DIFF
   INTEGER:: nitera, maxiter
   EXTERNAL:: fun

   OPEN(34, file='P3-23-24-res1.dat')
   preci=0.0001

   A=0.d0
   B=2.*4.*ATAN(1.d0)

   maxiter = NINT(LOG((B-A)/preci)/LOG(2.))+1
   WRITE(34, "(A35,I5)") "El nombre màxim d'iteracions és" , maxiter



   DO  nitera=1, maxiter

      valorarrel=(A+B)/2

      CALL fun(A, FA, valordf)
      CALL fun(B, FB, valordf)
      CALL fun(valorarrel, FC, valordf)

      IF (FA*FB .ge. 0.) THEN
         WRITE(34, "(A35,f20.12,f20.12,f20.12,f20.12)") 'La funció no canvia de signe.', A, B, FA, FB
         EXIT
      END IF

      IF (FC .eq. 0.) THEN
         WRITE(34, "(A35,f20.12)") 'Solució exacta: x=', valorarrel
         EXIT
      END IF
      IF (FA*FC .lt. 0.) THEN
         B=valorarrel
      ELSE
         A=valorarrel
      END IF
      DIFF=(B-A)

      IF (DIFF .le. PRECI) THEN
         WRITE(34, "(A35,f20.12)") 'Solució aproximada x=', valorarrel
         WRITE(34, "(A35,f20.12)") 'ERROR<', DIFF
         EXIT
      END IF
      WRITE(34, "(A35,f20.12,A35,f20.12,A35,f20.12)") 'Iteració número:', nitera, 'Valor arrel=', valorarrel, 'ERROR=', DIFF
   END DO

   CLOSE(34)
   RETURN
END
