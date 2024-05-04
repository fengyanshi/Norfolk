
!-------------------------------------------------------------------
!   This subroutine is used to pass coupling variables into ghost cells                                                         
!   Called by
!      TVD_SOLVER
!   Update: 04/28/2012 Fengyan Shi, University of Delaware   
!   1) used format I/O  2) read through read(11,119)                                      
!-------------------------------------------------------------------
SUBROUTINE OneWayCoupling
    USE GLOBAL
    IMPLICIT NONE

119      FORMAT(5E16.6)  ! this is a fixed format for I/O
 
! determine time slot

    IF(TIME>TIME_COUPLING_1.AND.TIME>TIME_COUPLING_2) THEN
         TIME_COUPLING_1=TIME_COUPLING_2

         READ(11,*,END=120) TIME_COUPLING_2 
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             READ(11,*,END=120)   ! east

             U_COUPLING_EAST(:,1)=U_COUPLING_EAST(:,2)
             V_COUPLING_EAST(:,1)=V_COUPLING_EAST(:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)

             READ(11,119,END=120)(U_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)(V_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
         ELSE
             READ(11,*,END=120)   ! east            
         ENDIF
! west
         IF(N_COUPLING_WEST.GT.0)THEN
             READ(11,*,END=120)   ! west

             U_COUPLING_WEST(:,1)=U_COUPLING_WEST(:,2)
             V_COUPLING_WEST(:,1)=V_COUPLING_WEST(:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)


             READ(11,119,END=120)(U_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)(V_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)

         ELSE
             READ(11,*,END=120)   ! west            
         ENDIF
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             READ(11,*,END=120)   ! south

             U_COUPLING_SOUTH(:,1)=U_COUPLING_SOUTH(:,2)
             V_COUPLING_SOUTH(:,1)=V_COUPLING_SOUTH(:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)


             READ(11,119,END=120)(U_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)(V_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
         ELSE
             READ(11,*,END=120)   ! south            
         ENDIF
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             READ(11,*,END=120)   ! north

             U_COUPLING_NORTH(:,1)=U_COUPLING_NORTH(:,2)
             V_COUPLING_NORTH(:,1)=V_COUPLING_NORTH(:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)

             READ(11,119,END=120)(U_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)(V_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
         ELSE
             READ(11,*,END=120)   ! north            
         ENDIF


    ENDIF  ! time>time_2 and time_1

120 CONTINUE

    tmp2=ZERO
    tmp1=ZERO

    IF(TIME>TIME_COUPLING_1)THEN
      IF(TIME_COUPLING_1.EQ.TIME_COUPLING_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=ZERO
      ELSE
      tmp2=(TIME_COUPLING_2-TIME) &
            /MAX(SMALL, ABS(TIME_COUPLING_2-TIME_COUPLING_1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1


! west boundary
   IF(N_COUPLING_WEST>0)THEN



     IF(IN_DOMAIN_WEST)THEN
      DO J=Kstart_WEST,Kend_WEST 
      DO I=1,Nghost
        ETA(I,J)=Z_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +Z_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2




        U(I,J)=U_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +U_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        V(I,J)=V_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +V_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2

        HU(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))
        HV(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))
        Ubar(I,J)=Jaco(I,J)*U(I,J)*(Depth(I,J)+ETA(I,J))
        Vbar(I,J)=Jaco(I,J)*V(I,J)*(Depth(I,J)+ETA(I,J))
      ENDDO
      ENDDO
     ENDIF  ! end in domain



    ENDIF ! end of n_coupling_west>0

! east boundary
   IF(N_COUPLING_EAST>0)THEN



     IF(IN_DOMAIN_EAST)THEN
      DO J=Kstart_EAST,Kend_EAST  
      DO I=Iend+1,Iend+Nghost
        ETA(I,J)=Z_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +Z_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2




        U(I,J)=U_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +U_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        V(I,J)=V_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +V_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2


        HU(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))
        HV(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))
        Ubar(I,J)=Jaco(I,J)*U(I,J)*(Depth(I,J)+ETA(I,J))
        Vbar(I,J)=Jaco(I,J)*V(I,J)*(Depth(I,J)+ETA(I,J))
      ENDDO
      ENDDO
     ENDIF  ! end in domain



    ENDIF ! end of n_coupling_east>0

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN



     IF(IN_DOMAIN_SOUTH)THEN
      DO I=Kstart_SOUTH,Kend_SOUTH  
      DO J=1,Nghost
        ETA(I,J)=Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2




        U(I,J)=U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        V(I,J)=V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2

        HU(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))
        HV(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))
        Ubar(I,J)=Jaco(I,J)*U(I,J)*(Depth(I,J)+ETA(I,J))
        Vbar(I,J)=Jaco(I,J)*V(I,J)*(Depth(I,J)+ETA(I,J))
      ENDDO
      ENDDO
     ENDIF  ! end in domain



    ENDIF ! end of n_coupling_south>0

! north boundary
   IF(N_COUPLING_NORTH>0)THEN



     IF(IN_DOMAIN_NORTH)THEN
      DO I=Kstart_NORTH,Kend_NORTH  
      DO J=Jend+1,Jend+Nghost
        ETA(I,J)=Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2




        U(I,J)=U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        V(I,J)=V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2

        HU(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))
        HV(I,J)=(Depth(I,J)+ETA(I,J))*(U(I,J)*L12(I,J)+V(I,J)*L22(I,J))
        Ubar(I,J)=Jaco(I,J)*U(I,J)*(Depth(I,J)+ETA(I,J))
        Vbar(I,J)=Jaco(I,J)*V(I,J)*(Depth(I,J)+ETA(I,J))
      ENDDO
      ENDDO
     ENDIF  ! end in domain



    ENDIF ! end of n_coupling_north>0

END SUBROUTINE OneWayCoupling

! end coupling


!-------------------------------------------------------------------
!   This subroutine is used for open boundary conditions                                                        
!   Called by
!      MAIN
!   Update: 05/17/2011 Fengyan Shi, University of Delaware                                        
!-------------------------------------------------------------------
SUBROUTINE TIDE
   USE GLOBAL
    IMPLICIT NONE
    INTEGER :: Kpoi,Kcon,Ip,Jp

    IF(ETA_CLAMPED)THEN
      tmp2=TideStartDate*24.0_SP*3600.0_SP+TIME  ! real time
      DO Kpoi=1,NumEtaPoint
        tmp1=ZERO
        DO Kcon = 1,NumConstituent
          tmp1=tmp1+TideFac(Kcon)*TideAmp(Kcon,Kpoi)*COS(2.0_SP*pi/(TidePeriod(Kcon,Kpoi)*3600.0_SP)*tmp2 &
                -TidePha(Kcon,Kpoi)*DEG2RAD+Tideu0(Kcon)*DEG2RAD)
        ENDDO

        ETA(I_eta_bd(kpoi)+Nghost,J_eta_bd(Kpoi)+Nghost) = tmp1&
                *TANH(TIME/3600.0_SP/6.0_SP) 

      ENDDO

    ENDIF

END SUBROUTINE TIDE


!-------------------------------------------------------------------
!   This subroutine is used to collect data into ghost cells                                                         
!   Called by
!      TVD_SOLVER
!
!   Call PHI_COLL
!
!   Update: 07/09/2010 Fengyan Shi, 
!     1) use dummy variables 2) add vtype=3                                      
!-------------------------------------------------------------------
SUBROUTINE EXCHANGE
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: VTYPE
    REAL(SP),DIMENSION(Mloc,Nloc) :: rMASK

    VTYPE=1
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta,VTYPE,PERIODIC_X,PERIODIC_Y)

    rMASK = MASK ! for periodic boundary condition
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,rMASK,VTYPE,PERIODIC_X,PERIODIC_Y)  
    MASK = rMASK  

    VTYPE=2
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ubar,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Usurf,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ubott,VTYPE,PERIODIC_X,PERIODIC_Y)
    VTYPE=3
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vbar,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vsurf,VTYPE,PERIODIC_X,PERIODIC_Y)
    CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vbott,VTYPE,PERIODIC_X,PERIODIC_Y)
! etaR x mask is a wrong idea
!    Eta=Eta*MASK

    U=U*MASK
    V=V*MASK
    HU=HU*MASK
    HV=HV*MASK
    Ubar=Ubar*MASK
    Vbar=Vbar*MASK
    
END SUBROUTINE EXCHANGE

!-------------------------------------------------------------------
!   This subroutine is used to collect data into ghost cells
!   Called by
!      Exchange
!
!   Update: 09/07/2010 Fengyan Shi, fix:
!   1) u v symmetric problem, 2) remove use global 3) fix bug
!-------------------------------------------------------------------
SUBROUTINE PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,PERIODIC_X,PERIODIC_Y)
    USE PARAM



    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)
    LOGICAL, INTENT(IN) :: PERIODIC_X,PERIODIC_Y

! for Eta
    IF(VTYPE==1) THEN  ! for eta
      ! x-direction



     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(K,J)=PHI(Iend-Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
     ENDIF







     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend  
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
     ENDIF





      ! y-direction and corners



      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jend-Nghost+K)
        ENDDO
        ENDDO 
      ELSE   
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jbeg+Nghost-K)
        ENDDO
        ENDDO
      ENDIF







      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO 
      ELSE
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=PHI(I,Jend-K+1)
        ENDDO
        ENDDO
      ENDIF





! for u
    ELSEIF(VTYPE==2) THEN  ! for u (x-mirror condition)
      ! x-direction



     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(K,J)=PHI(Iend-Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
     ENDIF







     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO
     ENDIF





      ! y-direction and corners




      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jend-Nghost+K)
        ENDDO
        ENDDO 
      ELSE   
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jbeg+Nghost-K)
        ENDDO
        ENDDO
      ENDIF







      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO 
      ELSE
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=PHI(I,Jend-K+1)
        ENDDO
        ENDDO
      ENDIF





    ELSEIF(VTYPE==3) THEN ! for v (y-mirror condition)
! for v
      ! x-direction



     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(K,J)=PHI(Iend-Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO
     ENDIF







     IF(PERIODIC_X)THEN
      DO J=1,Nloc
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Nghost+K,J)
      ENDDO
      ENDDO
     ELSE
      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO
     ENDIF




      ! y-direction and corners




      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=PHI(I,Jend-Nghost+K)
        ENDDO
        ENDDO 
      ELSE   
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
        ENDDO
        ENDDO
      ENDIF







      IF(PERIODIC_Y)THEN
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=PHI(I,Nghost+K)
        ENDDO
        ENDDO 
      ELSE
        DO I=1,Mloc
        DO K=1,Nghost
          PHI(I,Jend+K)=-PHI(I,Jend-K+1)
        ENDDO
        ENDDO
      ENDIF




! for cross-derivatives
    ELSEIF(VTYPE==4) THEN ! VTYPE==4 for u and v cross-mirror
     ! x-direction



      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=ZERO
      ENDDO
      ENDDO







      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=ZERO
      ENDDO
      ENDDO




      ! y-direction and corners, this one is not an exact solution



      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=ZERO
      ENDDO
      ENDDO







      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=ZERO
      ENDDO
      ENDDO





! for symmetric
    ELSEIF(VTYPE==5)THEN
      ! x-direction



      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=PHI(Ibeg+Nghost-K,J)
       ENDDO
      ENDDO







      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=PHI(Iend-K+1,J)
      ENDDO
      ENDDO




      ! y-direction and corners




      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO







      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=PHI(I,Jend-K+1)
      ENDDO
      ENDDO





! for anti-symmetric
      ELSE
      ! x-direction



      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
      ENDDO
      ENDDO 







      DO J=Jbeg,Jend
      DO K=1,Nghost
        PHI(Iend+K,J)=-PHI(Iend-K+1,J)
      ENDDO
      ENDDO 




      ! y-direction and corners



      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
      ENDDO
      ENDDO     







      DO I=1,Mloc
      DO K=1,Nghost
        PHI(I,Jend+K)=-PHI(I,Jend-K+1)
      ENDDO
      ENDDO     





    ENDIF









END SUBROUTINE PHI_COLL

! ---------------------------------------------------
!    This is subroutine to provide boundary conditions at edges of domain
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE BOUNDARY_CONDITION
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Xi,Deps,utmp,vtmp,dataU,dataAngle
     INTEGER :: Ip,Jp,Kpoi,Kcon


! four sides of computational domain
    IF(PERIODIC_X)THEN
!    do nothing
    ELSE









!  coupling
   IF(IN_DOMAIN_WEST)THEN
     DO J=Jbeg,Kstart_WEST-1
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ibeg,J)
      Gx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ibeg,J)
      ENDDO

     DO J=Kend_WEST+1,Jend
      P(Ibeg,J)=ZERO
      Xi=EtaRxR(Ibeg,J)
      Deps=Depthx(Ibeg,J)
      Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ibeg,J)
      Gx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ibeg,J)
      ENDDO

    ENDIF ! end in_domain

! not coupling






!  coupling
   IF(IN_DOMAIN_EAST)THEN
     DO J=Jbeg,Kstart_EAST-1
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Iend1,J)
      Gx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Iend1,J)
     ENDDO
     DO J=Kend_EAST+1,Jend
      P(Iend1,J)=ZERO
      Xi=EtaRxL(Iend1,J)
      Deps=Depthx(Iend1,J)
      Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Iend1,J)
      Gx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Iend1,J)
     ENDDO
   ENDIF ! end in_domain
! not coupling










   ENDIF ! end periodic

! y direction
   IF(PERIODIC_Y)THEN
!   do nothing
   ELSE










   IF(IN_DOMAIN_SOUTH)THEN
     DO I=Ibeg,Kstart_SOUTH-1
      Q(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Fy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(I,Jbeg)
      ENDDO
     DO I=Kend_SOUTH+1,Iend
      Q(I,Jbeg)=ZERO
      Xi=EtaRyR(I,Jbeg)
      Deps=Depthy(I,Jbeg)
      Fy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(I,Jbeg)
      Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(I,Jbeg)
      ENDDO
   ENDIF ! end in_domain
! not coupling





   IF(IN_DOMAIN_NORTH)THEN
     DO I=Ibeg,Kstart_NORTH-1
      Q(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Fy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(I,Jend1)
     ENDDO
     DO I=Kend_NORTH+1,Iend
      Q(I,Jend1)=ZERO
      Xi=EtaRyL(I,Jend1)
      Deps=Depthy(I,Jend1)
      Fy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(I,Jend1)
      Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(I,Jend1)
     ENDDO
   ENDIF ! end in_domain
! not coupling










    ENDIF  ! end periodic

! mask points
! Jeff pointed out the loop should be Jbeg-1, Jend+1
! The problem is that the fluxes on the inter-processor boundaries may be
!modified if the point next to the boundary (e.g., in the ghost cells,
!managed by a different processor) is land, but as is the routine doesn't
!check for this. 

     DO j=Jbeg-1,Jend+1
     DO i=Ibeg-1,Iend+1
      IF(MASK(I,J)<1)THEN
        P(I,J)=ZERO
! Jeff reported a bug here for parallel version







        IF(I/=Ibeg)THEN

!new splitting method
      Xi=EtaRxL(I,J)
      Deps=Depthx(I,J)

         Fx(I,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I-1,J)*JL11x(I,J)
         Gx(I,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I-1,J)*JL12x(I,J)




        ELSE
         Fx(I,J)=ZERO
         Gx(I,J)=ZERO
        ENDIF

        P(I+1,J)=ZERO
! EAST

! Jeff also here







        IF(I/=Iend)THEN

! new splitting method
      Xi=EtaRxR(I+1,J)
      Deps=Depthx(I+1,J)

         Fx(I+1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I+1,J)*JL11x(I+1,J)
         Gx(I+1,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I+1,J)*JL12x(I+1,J)




        ELSE
         Fx(I+1,J)=ZERO
         Gx(I+1,J)=ZERO
        ENDIF


        Q(I,J)=ZERO
        Fy(I,J)=ZERO

! SOUTH
! Jeff also here







        IF(J/=Jbeg)THEN

! new splitting method
      Xi=EtaRyL(I,J)
      Deps=Depthy(I,J)

         Gy(I,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I,J-1)*JL22y(I,J)
         Fy(I,J)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I,J-1)*JL21y(I,J)




        ELSE
         Gy(I,J)=ZERO
         Fy(I,J)=ZERO
        ENDIF

        Q(I,J+1)=ZERO

! NORTH

! Jeff also here







        IF(J/=Jend)THEN

! new splitting method
      Xi=EtaRyR(I,J+1)
      Deps=Depthy(I,J+1)

         Gy(I,J+1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I,J+1)*JL22y(I,J+1)
         Fy(I,J+1)=0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*MASK(I,J+1)*JL21y(I,J+1)




        ELSE
         Gy(I,J+1)=ZERO
         FY(I,J+1)=ZERO
        ENDIF
      ENDIF
     ENDDO
     ENDDO

! FLUX --- USUALLY FOR RIVER DISCHARGE CONDITION

    IF(FLUX_CLAMPED)THEN
      IF(TIME>TimeFlux(icount_riverdata).AND.icount_riverdata<NumTimeData)THEN
        icount_riverdata=icount_riverdata+1
      ENDIF
      IF(icount_riverdata>1)THEN  ! river start discharge
       IF(TIME>TimeFlux(icount_riverdata))THEN
        tmp2=ZERO
       ELSE
        tmp2=(TimeFlux(icount_riverdata)-TIME) &
           /(TimeFlux(icount_riverdata)-TimeFlux(icount_riverdata-1))
       ENDIF
        DO Kpoi=1,NumFluxPoint




         Ip=I_flux_bd(Kpoi)+Nghost
         Jp=J_flux_bd(Kpoi)+Nghost






         IF(RIVER_ORIENTATION(Kpoi)=='W')THEN  !  left boundary

           dataU = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
           dataAngle=FluxAngle(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxAngle(Kpoi,icount_riverdata-1)*tmp2

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

           Xi=EtaRxR(Ip,Jp)
           Deps=Depthx(Ip,Jp)
! NOTE there was a big mistake here, P and Q should be contravariant 
! component, I corrected on 01/20/2012
!      the previous is     P(Ip,Jp) = utmp*JacoX(Ip,Jp)
           
           P(Ip,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp)

           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)
           
           Fx(Ip,Jp) = P(Ip,Jp)*utmp  &
                  +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ip,Jp)

           Gx(Ip,Jp) = P(Ip,Jp)*vtmp  &
                  +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ip,Jp)



         ELSEIF(RIVER_ORIENTATION(Kpoi)=='E')THEN  ! right boundary

           dataU = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
           dataAngle=FluxAngle(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxAngle(Kpoi,icount_riverdata-1)*tmp2

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

           Xi=EtaRxL(Ip+1,Jp)
           Deps=Depthx(Ip+1,Jp) 


           P(Ip+1,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip+1,Jp)

!           P(Ip+1,Jp)=utmp*JacoX(Ip+1,Jp)

           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)

           Fx(Ip+1,Jp) = P(Ip+1,Jp)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ip+1,Jp)
           Gx(Ip+1,Jp) = P(Ip+1,Jp)*vtmp  &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ip+1,Jp)









         ELSE
           ! do nothing for x direction   
         ENDIF

         IF(RIVER_ORIENTATION(Kpoi)=='S')THEN  !left boundary

           dataU = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
           dataAngle=FluxAngle(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxAngle(Kpoi,icount_riverdata-1)*tmp2

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

!           utmp = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
!                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
 
           Xi=EtaRyR(Ip,Jp)
           Deps=Depthy(Ip,Jp)

           P(Ip,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp)
           Q(Ip,Jp) = (utmp*L12(Ip,Jp)+vtmp*L22(Ip,Jp))*JacoY(Ip,Jp)

! change          Q(Ip,Jp)=utmp*JacoY(Ip,Jp)


           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)

           Gy(Ip,Jp) = Q(Ip,Jp)*vtmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(Ip,Jp)
           Fy(Ip,Jp) = Q(Ip,Jp)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(Ip,Jp)









          ELSEIF(RIVER_ORIENTATION(Kpoi)=='N')THEN ! right boundary

           dataU = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
           dataAngle=FluxAngle(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
                       +FluxAngle(Kpoi,icount_riverdata-1)*tmp2

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)


!           utmp = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
!                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2

           Xi=EtaRyL(Ip,Jp+1)
           Deps=Depthy(Ip,Jp+1)  

           P(Ip,Jp+1) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp+1)
           Q(Ip,Jp+1) = (utmp*L12(Ip,Jp)+vtmp*L22(Ip,Jp))*JacoY(Ip,Jp+1)


!           Q(Ip,Jp+1)=utmp

           utmp=utmp/MAX(Deps,SMALL)  
           vtmp=vtmp/MAX(Deps,SMALL)
              
           Gy(Ip,Jp+1) = Q(Ip,Jp+1)*vtmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(Ip,Jp+1) 
           Fy(Ip,Jp+1) = Q(Ip,Jp+1)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(Ip,Jp+1)









          ELSE
            ! do nothing for y direction
          ENDIF




        ENDDO ! end Kpoi

      ENDIF     ! end river start discharge

    ELSEIF(FLUX_TIDE)THEN  ! end flux_clamped and start tidal flux boundary

! !!!!!!!!!!!!!
        tmp2=TideStartDate*24.0_SP*3600.0_SP+TIME  ! real time
        DO Kpoi=1,NumFluxPoint




         Ip=I_flux_bd(Kpoi)+Nghost
         Jp=J_flux_bd(Kpoi)+Nghost






         IF(FLUX_ORIENTATION(Kpoi)=='W')THEN  !  left boundary

           dataU=ZERO
           DO Kcon=1,NumConstituent
             dataU = dataU + FluxU_tide(Kcon,Kpoi)*COS(2.0_SP*pi/(FluxPeriod(Kcon,Kpoi)*3600.0_SP)*tmp2 &
                -FluxPha(Kcon,Kpoi)*DEG2RAD)
           ENDDO

           dataAngle=FluxAngle_tide(Kpoi)  !$$$$

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

           Xi=EtaRxR(Ip,Jp)
           Deps=Depthx(Ip,Jp)
           
           P(Ip,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp)

           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)
           
           Fx(Ip,Jp) = P(Ip,Jp)*utmp  &
                  +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ip,Jp)

           Gx(Ip,Jp) = P(Ip,Jp)*vtmp  &
                  +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ip,Jp)



         ELSEIF(FLUX_ORIENTATION(Kpoi)=='E')THEN  ! right boundary


           dataU=ZERO
           DO Kcon=1,NumConstituent
             dataU = dataU + FluxU_tide(Kcon,Kpoi)*COS(2.0_SP*pi/(FluxPeriod(Kcon,Kpoi)*3600.0_SP)*tmp2 &
                -FluxPha(Kcon,Kpoi)*DEG2RAD)
           ENDDO

           dataAngle=FluxAngle_tide(Kpoi)

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

           Xi=EtaRxL(Ip+1,Jp)
           Deps=Depthx(Ip+1,Jp) 


           P(Ip+1,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip+1,Jp)

           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)

           Fx(Ip+1,Jp) = P(Ip+1,Jp)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL11x(Ip+1,Jp)
           Gx(Ip+1,Jp) = P(Ip+1,Jp)*vtmp  &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL12x(Ip+1,Jp)


         ELSE
           ! do nothing for x direction   
         ENDIF

         IF(FLUX_ORIENTATION(Kpoi)=='S')THEN  !left boundary

           dataU=ZERO
           DO Kcon=1,NumConstituent
             dataU = dataU + FluxU_tide(Kcon,Kpoi)*COS(2.0_SP*pi/(FluxPeriod(Kcon,Kpoi)*3600.0_SP)*tmp2 &
                -FluxPha(Kcon,Kpoi)*DEG2RAD)
           ENDDO

           dataAngle=FluxAngle_tide(Kpoi)

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)

!           utmp = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
!                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2
 
           Xi=EtaRyR(Ip,Jp)
           Deps=Depthy(Ip,Jp)

           P(Ip,Jp) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp)
           Q(Ip,Jp) = (utmp*L12(Ip,Jp)+vtmp*L22(Ip,Jp))*JacoY(Ip,Jp)

! change          Q(Ip,Jp)=utmp*JacoY(Ip,Jp)


           utmp=utmp/MAX(Deps,SMALL)
           vtmp=vtmp/MAX(Deps,SMALL)

           Gy(Ip,Jp) = Q(Ip,Jp)*vtmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(Ip,Jp)
           Fy(Ip,Jp) = Q(Ip,Jp)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(Ip,Jp)


          ELSEIF(FLUX_ORIENTATION(Kpoi)=='N')THEN ! right boundary

           dataU=ZERO
           DO Kcon=1,NumConstituent
             dataU = dataU + FluxU_tide(Kcon,Kpoi)*COS(2.0_SP*pi/(FluxPeriod(Kcon,Kpoi)*3600.0_SP)*tmp2 &
                -FluxPha(Kcon,Kpoi)*DEG2RAD)
           ENDDO

           dataAngle=FluxAngle_tide(Kpoi)

           utmp=dataU*COS(dataAngle*Deg2Rad)
           vtmp=dataU*SIN(dataAngle*Deg2Rad)


!           utmp = FluxU_bd(Kpoi,icount_riverdata)*(1.0_SP-tmp2)&
!                       +FluxU_bd(Kpoi,icount_riverdata-1)*tmp2

           Xi=EtaRyL(Ip,Jp+1)
           Deps=Depthy(Ip,Jp+1)  

           P(Ip,Jp+1) = (utmp*L11(Ip,Jp)+vtmp*L12(Ip,Jp))*JacoX(Ip,Jp+1)
           Q(Ip,Jp+1) = (utmp*L12(Ip,Jp)+vtmp*L22(Ip,Jp))*JacoY(Ip,Jp+1)


!           Q(Ip,Jp+1)=utmp

           utmp=utmp/MAX(Deps,SMALL)  
           vtmp=vtmp/MAX(Deps,SMALL)
              
           Gy(Ip,Jp+1) = Q(Ip,Jp+1)*vtmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL22y(Ip,Jp+1) 
           Fy(Ip,Jp+1) = Q(Ip,Jp+1)*utmp &
                 +0.5_SP*GRAV*(Xi*Xi+2.0_SP*Xi*Deps)*JL21y(Ip,Jp+1)


          ELSE
            ! do nothing for y direction
          ENDIF




        ENDDO ! end Kpoi

! !!!!!!!!!!!!
    ENDIF ! end flux_tidal

END SUBROUTINE BOUNDARY_CONDITION


! --------------------------------------------------
!    This is subroutine to provide solitary wave at left boundary
!    it can be specified in input.txt giving 'SOL'
!    called by
!       - MAIN
!    Last Update: 05/28/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE SOLITARY_WAVE_LEFT_BOUNDARY
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP):: aa,h00,c1,tex,tlag,zz
     INTEGER::Iwavemaker

       Iwavemaker=Ibeg
       aa=AMP_SOLI
       h00=DEP_SOLI
       c1=sqrt(GRAV*h00*(1.0_SP+aa/h00))
       DO J=1,Nloc
         tex=sqrt(0.75_SP*aa/h00**3)
         tlag=4.0_SP*h00/sqrt(aa/h00)
         zz=aa/COSH(tex*(Lag_soli-c1*TIME))**2
         Eta(Iwavemaker,J)=zz
         H(Iwavemaker,J)=Eta(Iwavemaker,J)+Depth(Iwavemaker,J) 
! note: can not provide u and hu at boundary for dispersive equations!
!         U(Iwavemaker,J)= SQRT(grav/h00)*zz
!         HU(Iwavemaker,J)=h00*U(Iwavemaker,J)       
       enddo   
     
END SUBROUTINE SOLITARY_WAVE_LEFT_BOUNDARY
