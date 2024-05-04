
! ---------------------------------------------------
!    This is subroutine to solve PERIODIC tridiagonal matrix
!    Last Update: 01/17/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
!solve the following cyclic tridiagonal equation
!*************************************************************!
!*                                                           *!
!*         (B1 C1                      A1 )   (x1)     (b1)  *!
!*         (A2 B2 C2                      )   (x2)     (b2)  *!
!*         (   A3 B3 C3                   )   (x3)     (b3)  *!
!*         (      A4 B4 C4                )   (x4)====       *!
!*         (         A5 B5  C5            )   ... ====  ...  *!
!*         (            ... ... ...       )   ...       ...  *!
!*         (                AN-1 BN-1 CN-1)   (xN-1)   (bN-1)*!
!*         (CN                   AN   BN  )   (xN)     (bN)  *!
!*                                                           *!
!*                                                           *!
!*************************************************************!
!    where A are alpha, B are beta, C are gama
!  CALL TRIG_PERIODIC(Aper,Bper,Cper,Dper,Vper,Jend-Jbeg+1)

SUBROUTINE TRIG_PERIODIC (alpha, beta,gama, b, x, N)  
         USE PARAM
         IMPLICIT NONE
!----dummy variables--
         INTEGER,INTENT(IN) :: N
         REAL(SP),DIMENSION(N),INTENT(IN) :: alpha, beta,gama,b
         REAL(SP),DIMENSION(N),INTENT(OUT) ::  x
!----local variables--
          REAL(SP),DIMENSION(N) :: betaPrime,u,z
          REAL(SP) :: c,fact
          INTEGER :: II

          c=-beta(1)   !the minus sign makes sure betaPrime(1) non zero

          betaPrime(1) = beta(1) - c     !the first tirdiagonal element
          do II=2,N-1,1
           betaPrime(II)=beta(II)
          enddo
   
          betaPrime(N) = beta(N) - alpha(1) * gama(N) / c    !the last tridiagonal element

          call tri_ge(alpha, betaPrime, gama, b, x, N)       !solve for x  !first equation (13) on notes page 148-4-3

          u(1) = c;                                          !first u,
          do II=2,N-1
           u(II)=0.d0
          enddo
          u(N) = gama(N)                                    !last U is same as alpha

        call tri_ge(alpha, betaPrime, gama, u, z, N)       !solve for z  !second equation (13) on notes page 148-4-3

        fact = (x(1) + alpha(1) * x(N) / c)   &
             / (1.0 + z(1) + alpha(1) * z(N) / c)

          do II=1,N
             x(II) =x(II)- fact * z(II)       !construct final results
          enddo

END SUBROUTINE TRIG_PERIODIC




SUBROUTINE TRI_GE(alpha,beta,gama,b, x, N)
!----basically same as subroutine trig()----but allows diagonal variables not equal to unit
!*************************************************************!
!*                                                           *!
!*         (B1 C1                         )   (x1)     (b1)  *!
!*         (A2 B2 C2                      )   (x2)     (b2)  *!
!*         (   A3 B3 C3                   )   (x3)     (b3)  *!
!*         (      A4 B4 C4                )   (x4)====       *!
!*         (         A5 B5  C5            )   ... ====  ...  *!
!*         (            ... ... ...       )   ...       ...  *!
!*         (                An-1 Bn-1 Cn-1)   (xn-1)   (bn-1)*!
!*         (                     An   Bn  )   (xn)     (bn)  *!
!*                                                           *!
!*                                                           *!
!*************************************************************!
! where A are alpha, B are beta, C are gama
        USE PARAM
        IMPLICIT NONE
!----dummy variables---
        INTEGER,INTENT(IN) :: N
        REAL(SP),DIMENSION(N),INTENT(IN) :: alpha,beta, gama, b
        REAL(SP),DIMENSION(N),INTENT(OUT) :: x
!-----local variables
        REAL(SP),DIMENSION(N) :: betaPrime, bPrime
        REAL(SP) :: coeff
        INTEGER :: II
!Perform forward elimination
         betaPrime(1) = beta(1)
         bPrime(1) = b(1)
 
         do II=2,N
               coeff = alpha(II) / betaPrime(II-1)
               betaPrime(II) = beta(II) - coeff * gama(II-1)
               bPrime(II) = b(II) - coeff * bPrime(II-1)
         enddo

! Perform back substitution
         x(N) = bPrime(N) / betaPrime(N)
         do II=N-1,1,-1
              x(II) = (bPrime(II) - gama(II) * x(II+1)) / betaPrime(II)
         enddo

END SUBROUTINE TRI_GE



! ---------------------------------------------------
!    This is subroutine to solve tridiagonal matrix
!    Last Update: 09/24/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE TRIG(A,C,D,Z,M,Ibeg,Iend)
        USE PARAM
        IMPLICIT NONE
        INTEGER :: II
        INTEGER, INTENT(IN) :: M,Ibeg,Iend
        REAL(SP),DIMENSION(M),INTENT(INOUT) :: A,C,D
        REAL(SP),DIMENSION(M),INTENT(OUT) :: Z

        DO II=Ibeg+1,Iend
          IF(A(II).NE.ZERO)THEN
            C(II)=C(II)/A(II)/(1.0_SP/A(II)-C(II-1))
            D(II)=(D(II)/A(II)-D(II-1))/(1.0_SP/A(II)-C(II-1))
          ENDIF
        ENDDO

        Z(Iend)=D(Iend)

        DO II=Iend-1,Ibeg,-1
          Z(II)=D(II)-C(II)*Z(II+1)
        ENDDO

END SUBROUTINE TRIG














