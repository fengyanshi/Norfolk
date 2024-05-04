

! --------------------------------------------------
!    This is subroutine to solve 3d profile
!    called by
!         MAIN
!    Last Update: 08/81/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE SOLVE_3D
    USE PARAM
    USE GLOBAL, ONLY : U,V,Mloc,Nloc,Ibeg,Iend,Jbeg,Jend, &
                       nu_total,H,Tau_bx,Tau_by,Tau_sx,Tau_sy,&
                       WindForce,dd11,ee11,ff11,dd12,ee12,ff12, &
                       Ubott,Vbott,MinDepthFrc,D_11,D_12,D_22, &
                       tmp4preview,DISP_EMPIRICAL,Usurf,Vsurf
    USE PASS, ONLY :  WaveFxSC,WaveFySC,WaveFluxXSC,WaveFluxYSC, &
                      WaveAngleSC,WaveHeightSC
    IMPLICIT NONE
    REAL(SP) :: AAx,AAy,rate

! preprocess for case with less components
    IF(.NOT.ALLOCATED(dd11))THEN
      ALLOCATE(dd11(Mloc,Nloc),dd12(Mloc,Nloc),ee11(Mloc,Nloc),ee12(Mloc,Nloc),&
               ff11(Mloc,Nloc),ff12(Mloc,Nloc) &
              )
    ENDIF


    IF(.NOT.WindForce)THEN
      Tau_sx=ZERO
      Tau_sy=ZERO
    ENDIF
    IF(.NOT.ALLOCATED(WaveFxSC))THEN
       ALLOCATE(WaveFxSC(Mloc,Nloc),WaveFySC(Mloc,Nloc))
       WaveFxSC=ZERO
       WaveFySC=ZERO
    ENDIF
    IF(.NOT.ALLOCATED(WaveFluxXSC))THEN
       ALLOCATE(WaveFluxXSC(Mloc,Nloc),WaveFluxYSC(Mloc,Nloc))
       WaveFluxXSC=ZERO
       WaveFluxYSC=ZERO
    ENDIF
    IF(.NOT.ALLOCATED(WaveAngleSC))THEN
       ALLOCATE(WaveAngleSC(Mloc,Nloc))
       WaveAngleSC=ZERO
    ENDIF
    IF(.NOT.ALLOCATED(WaveHeightSC))THEN
           ALLOCATE(WaveHeightSC(Mloc,Nloc))
           WaveHeightSC = ZERO
    ENDIF



     DO J=Jbeg,Jend
       DO I=Ibeg,Iend
         AAx=(Tau_bx(I,J)-Tau_sx(I,J))/MAX(nu_total(I,J),SMALL)
         AAy=(Tau_by(I,J)-Tau_sy(I,J))/MAX(nu_total(I,J),SMALL)

         ff11(I,J)=-1.0_SP/3.0_SP*AAx*H(I,J) &
              - WaveFluxXSC(I,J)/MAX(H(I,J),MinDepthFrc)
         ff12(I,J)=-1.0_SP/3.0_SP*AAy*H(I,J) &
              - WaveFluxYSC(I,J)/MAX(H(I,J),MinDepthFrc)

         ee11(I,J)=AAx
         ee12(I,J)=AAy

         dd11(I,J)=-AAx*0.5_SP/MAX(H(I,J),MinDepthFrc)
         dd12(I,J)=-AAy*0.5_SP/MAX(H(I,J),MinDepthFrc)

         Ubott(I,J)=ff11(I,J) + U(I,J)
         Vbott(I,J)=ff12(I,J) + V(I,J)

         Usurf(I,J)=ff11(I,J)+dd11(I,J)*H(I,J)*H(I,J)+ &
                    ee11(I,J)*H(I,J)+U(I,J) + &
                   WaveFluxXSC(I,J)/MAX(H(I,J),MinDepthFrc)
         Vsurf(I,J)=ff12(I,J)+dd12(I,J)*H(I,J)*H(I,J)+ &
                    ee12(I,J)*H(I,J)+V(I,J) +&
                   WaveFluxYSC(I,J)/MAX(H(I,J),MinDepthFrc)

         IF(DISP_EMPIRICAL)THEN
! based on Putrevu and Svendsen's empirical formula 
           rate=WaveHeightSC(I,J)/Max(H(I,J),MinDepthFrc)
           D_11(I,J)=rate*rate
           IF(D_11(I,J).GT.0.5_SP)D_11(I,J)=0.5_SP
           D_11(I,J)=0.60_SP*D_11(I,J)
           D_22(I,J)=D_11(I,J)
         ELSE
! based on Putrevu and Svendsen's theoretical formula
! dispersion terms (these terms are greatly simplified, see notes)
         IF(H(I,J).GT.ZERO)THEN
         D_11(I,J)=1.0_SP/Max(nu_total(I,J),SMALL)*( &
            1.0_SP/63.0_SP*dd11(I,J)*dd11(I,J)*(H(I,J))**6  & 
           +1.0_SP/36.0_SP*(dd11(I,J)*ee11(I,J)+dd11(I,J)*ee11(I,J))  &
           *(H(I,J))**5  &
           +(1.0_SP/15.0_SP*dd11(I,J)*ff11(I,J)+1.0_SP/15.0_SP*dd11(I,J)*ff11(I,J)  &
           +1.0_SP/20.0_SP*ee11(I,J)*ee11(I,J))  &
             *(H(I,J))**4  &
           +1.0_SP/8.0_SP*(ee11(I,J)*ff11(I,J)+ee11(I,J)*ff11(I,J)) &
             *(H(I,J))**3 &
           +1.0_SP/3.0_SP*ff11(I,J)*ff11(I,J)*(H(I,J))**2  &
              )
	D_22(I,J)=1.0_SP/Max(nu_total(i,j),SMALL)*( &
            1.0_SP/63.0_SP*dd12(I,J)*dd12(I,J)*(H(I,J))**6 & 
           +1.0_SP/36.0_SP*(dd12(I,J)*ee12(I,J)+dd12(I,J)*ee12(I,J)) &
             *(H(I,J))**5  &
           +(1.0_SP/15.0_SP*dd12(I,J)*ff12(I,J)+1.0_SP/15.0_SP*dd12(I,J)*ff12(I,J) &
           +1.0_SP/20.0_SP*ee12(I,J)*ee12(I,J))  &
             *(H(I,J))**4  &
           +1.0_SP/8.0_SP*(ee12(I,J)*ff12(I,J)+ee12(I,J)*ff12(I,J))  &
             *(H(I,J))**3  &
           +1.0_SP/3.0_SP*ff12(I,J)*ff12(I,J)*(H(I,J))**2  &
               )
         ELSE
         D_11(I,J)=ZERO
         D_22(I,J)=ZERO
         ENDIF
        ENDIF ! end empirical
       ENDDO
     ENDDO


END SUBROUTINE SOLVE_3D

! --------------------------------------------------
!    This is subroutine to update wind
!    called by
!         MAIN
!    Last Update: 08/81/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE WIND_SPEED
    USE GLOBAL
    IMPLICIT NONE
      IF(TIME>TimeWind(icount_winddata).AND. &
               icount_winddata<NumTimeWindData)THEN
        icount_winddata=icount_winddata+1
      ENDIF
      IF(icount_winddata>1)THEN  ! wind start
       IF(TIME>TimeWind(icount_winddata))THEN
        tmp2=ZERO
       ELSE
        tmp2=(TimeWind(icount_winddata)-TIME) &
           /(TimeWind(icount_winddata)-TimeWind(icount_winddata-1))
       ENDIF
       WindU2D=WU(icount_winddata)*(1.0_SP-tmp2)&
                       +WU(icount_winddata-1)*tmp2
       WindV2D=WV(icount_winddata)*(1.0_SP-tmp2)&
                       +WV(icount_winddata-1)*tmp2
      ENDIF
END SUBROUTINE WIND_SPEED

! --------------------------------------------------
!    This is subroutine to update mask
!    note that mask also be updated in fluxes routine
!    called by
!         MAIN
!    Last Update: 12/06/2013 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE UPDATE_MASK
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::left,right,top,bottom
     REAL(SP)::FloodLimit
     INTEGER :: MASK4

! for the serial code, MASK at ghost cells keep no change






! Jeff did the following loop, also work for serial
!     DO J=Jbeg,Jend
!     DO I=Ibeg,Iend

     DO J=Jbeg-2,Jend+2
     DO I=Ibeg-2,Iend+2
! flood
     IF(MASK_STRUC(I,J)==1)THEN
       IF(MASK(I,J)<1)THEN
        MASK4=0
        FloodLimit = -Depth(I,j)+MinDepth
         ! left
        IF(I/=1)THEN
         IF(MASK(I-1,J)==1.AND.Eta(I-1,J)>FloodLimit.AND.U(I-1,J)>ZERO)THEN
           MASK(I,J)=1
         ENDIF
        ENDIF
         ! right
        IF(I/=Mloc)THEN
         IF(MASK(I+1,J)==1.AND.Eta(I+1,J)>FloodLimit.AND.U(I+1,J)<ZERO)THEN
           MASK(I,J)=1
         ENDIF
        ENDIF
         ! bottom
        IF(J/=1)THEN
         IF(MASK(I,J-1)==1.AND.Eta(I,J-1)>FloodLimit.AND.V(I,J-1)>ZERO)THEN
           MASK(I,J)=1
         ENDIF
        ENDIF
         ! top
        IF(J/=Nloc)THEN
         IF(MASK(I,J+1)==1.AND.Eta(I,J+1)>FloodLimit.AND.V(I,J+1)<ZERO)THEN
           MASK(I,J)=1
         ENDIF
         
        ENDIF
100   continue
! drying
       ELSE
         IF(Eta(I,J)<-Depth(I,J))THEN
          MASK(I,J)=0
!          Eta(I,J)=-MinDepth-Depth(I,J)    
         ENDIF    
       ENDIF
      ENDIF

! to avoid extreme depth gradient caused by depthx and depthy which were not
! treated in initialization, I reset depthx and depthy when drying 
! 01/19/2012

        IF(MASK(I,J)<1)THEN
         DepthX(I,J)=Depth(I-1,J)
         DepthX(I+1,J)=Depth(I+1,J)
         DepthY(I,J)=Depth(I,J-1)
         DepthY(I,J+1)=Depth(I,J+1)
        ENDIF    

     ENDDO
     ENDDO

! Jeff also did this loop
!     DO J=Jbeg,Jend
!     DO I=Ibeg,Iend
     DO J=Jbeg-1,Jend+1
     DO I=Ibeg-1,Iend+1
      MASK9(I,J)=MASK(I,J)*MASK(I-1,J)*MASK(I+1,J)  &
                *MASK(I+1,J+1)*MASK(I,J+1)*MASK(I-1,J+1) &
                *MASK(I+1,J-1)*MASK(I,J-1)*MASK(I-1,J-1) 
      IF(ABS(Eta(I,J))/MAX(DEPTH(I,J),MinDepthFrc)>SWE_ETA_DEP)THEN
       MASK9(I,J)=ZERO
      ENDIF

     ENDDO
     ENDDO
  





END SUBROUTINE UPDATE_MASK

! --------------------------------------------------
!    This is subroutine for all source terms including slope term dispersion É
!    called by
!       MAIN
!    Last Update: 05/28/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE SourceTerms
     USE GLOBAL
     USE PASS
     IMPLICIT NONE
     REAL(SP) :: diffusion
! local variables
     REAL(SP) :: HL11,HL12,HL22, &
                 HL11_1,HL11_2,HL12_1,HL12_2,HL22_1,HL22_2, &
                 NJL11,NJL12,NJL22, &
                 u_1,u_2,v_1,v_2, &
                 u_11,u_12,u_22,v_11,v_12,v_22, &
                 e11,e12,e22, &
                 Uu0

! for shorecirc-only-run
     IF(.NOT.ALLOCATED(WaveUbottSC))THEN
      ALLOCATE(WaveUbottSC(Mloc,Nloc))
      WaveUbottSC=ZERO
     ENDIF
     IF(.NOT.ALLOCATED(WaveAngleSC))THEN
       ALLOCATE(WaveAngleSC(Mloc,Nloc))
       WaveAngleSC=ZERO
     ENDIF

! depth gradient term
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend

! second order, move the second term to left-hand side


! for mixing terms
       u_1=(U(I+1,J)-U(I-1,J))*0.5_SP
       u_2=(U(I,J+1)-U(I,J-1))*0.5_SP
       v_1=(V(I+1,J)-V(I-1,J))*0.5_SP
       v_2=(V(I,J+1)-V(I,J-1))*0.5_SP
       u_11=U(I+1,J)-2.0_SP*U(I,J)+U(I-1,J)
       v_11=V(I+1,J)-2.0_SP*V(I,J)+V(I-1,J)
       u_22=U(I,J+1)-2.0_SP*U(I,J)+U(I,J-1)
       v_22=V(I,J+1)-2.0_SP*V(I,J)+V(I,J-1)
       u_12=0.25_SP*(U(I+1,J+1)-U(I-1,J+1)-U(I+1,J-1)+U(I-1,J-1))
       v_12=0.25_SP*(V(I+1,J+1)-V(I-1,J+1)-V(I+1,J-1)+V(I-1,J-1))

! smagorinsky
       e11=L11(I,J)*u_1+L12(I,J)*u_2
       e22=L22(I,J)*v_2+L12(I,J)*v_1
       e12=0.5_SP*(L12(I,J)*u_1+L22(I,J)*u_2+L11(I,J)*v_1+L12(I,J)*v_2)
      ! use 2-dimension just for analysis
       nu_smg(I,J)=Jaco(I,J)*C_smg*C_smg  &
                    *SQRT(2.0_SP*(e11*e11+2.0_SP*e12*e12+e22*e22))
!       nu_total(I,J)=nu_0+nu_bkgd+nu_smg(I,J)
! something strange when using smagorinsky, a lot of noises
! I remove smag for now
       nu_total(I,J)=nu_0+nu_bkgd

       HL11=H(I,J)*L11(I,J)
       HL12=H(I,J)*L12(I,J)
       HL22=H(I,J)*L22(I,J)
       HL11_1=0.5_SP*(H(I+1,J)*L11(I+1,J)-H(I-1,J)*L11(I-1,J))
       HL12_1=0.5_SP*(H(I+1,J)*L12(I+1,J)-H(I-1,J)*L12(I-1,J))
       HL22_1=0.5_SP*(H(I+1,J)*L22(I+1,J)-H(I-1,J)*L22(I-1,J))
       HL11_2=0.5_SP*(H(I,J+1)*L11(I,J+1)-H(I,J-1)*L11(I,J-1))
       HL12_2=0.5_SP*(H(I,J+1)*L12(I,J+1)-H(I,J-1)*L12(I,J-1))
       HL22_2=0.5_SP*(H(I,J+1)*L22(I,J+1)-H(I,J-1)*L22(I,J-1))
       NJL11=Jaco(I,J)*L11(I,J)
       NJL12=Jaco(I,J)*L12(I,J)
       NJL22=Jaco(I,J)*L22(I,J)

       SourceX(I,J)=GRAV*(Eta(I,J))*  &
          ((Depthx(I+1,J)*JL11x(I+1,J)-Depthx(I,J)*JL11x(I,J))/DX &
          +(Depthy(I,J+1)*JL21y(I,J+1)-Depthy(I,J)*JL21y(I,J))/DY)* &
          MASK(I,J) 

! ANNA start









! ANNA end

! friction
! curv-fitting beta1 and beta2, only 1 data for Vb=u0 and nu=0
       beta1=1.0_SP
       beta2=0.5_SP
     IF(TIME.GE.WC_LAG)THEN 
       Uu0=MAX(SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J)),WaveUbottSC(I,J))
     ELSE
       Uu0=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
     ENDIF
! calculate Cd using Manning formula
!      Cd=g*Manning^2/H^(1/3)

       IF(FRC_MANNING_DATA)THEN

        Cd=2.0_SP*grav*CD_2D(I,J)**2/(Max(H(I,J), MinDepthFrc))**(1.0_SP/3.0_SP)

       ELSE IF(Manning.gt.ZERO)THEN
        Cd=2.0_SP*grav*Manning**2/(Max(H(I,J), MinDepthFrc))**(1.0_SP/3.0_SP)
       ENDIF
       
       Tau_bx(I,J)=0.5_SP*Cd*(beta1*U(I,J) &
           +beta2*WaveUbottSC(I,J)*COS(WaveAngleSC(I,J)*Deg2Rad))*Uu0

       SourceX(I,J)=SourceX(I,J) &
          -Jaco(I,J)*Tau_bx(I,J) &
                        ! Coriolis
                    +Coriolis(I,J)*Vbar(I,J) 

! ANNA start











! ANNA end

! radiation stresses 
       SourceX(I,J)=SourceX(I,J) &
          +Jaco(I,J)*grav*WaveFxSC(I,J)

! ANNA start









! ANNA end

! lateral mixing terms and 3d dispersion (Disp = 1.0 set in init.fsc)

       diffusion =-(nu_total(I,J)-Disp*D_11(I,J))  &
               *NJL11*(2.0_SP*HL11*u_11+2.0_SP*u_1*HL11_1  &
               +2.0_SP*HL12*u_12+2.0_SP*u_2*HL12_1 )  &
          -(nu_total(I,J)-Disp*D_11(I,J))  &
               *NJL12*(2.0_SP*HL11*u_12+2.0_SP*u_1*HL11_2  &
               +2.0_SP*HL12*u_22+2.0_SP*u_2*HL12_2 )  &
          -(nu_total(I,J)-Disp*D_22(I,J))  &
               *NJL12*(HL12*u_11+u_1*HL12_1+HL22*u_12+u_2*HL22_1)  &
          -(nu_total(I,J)-Disp*D_11(I,J))  &
               *NJL12*(HL11*v_11+v_1*HL11_1+HL12*v_12+v_2*HL12_1)  &
          -(nu_total(I,J)-Disp*D_22(I,J))  &
               *NJL22*(HL12*u_12+u_1*HL12_2+HL22*u_22+u_2*HL22_2)  &
          -(nu_total(I,J)-Disp*D_11(I,J))  &
               *NJL22*(HL11*v_12+v_1*HL11_2+HL12*v_22+v_2*HL12_2 )
        SourceX(I,J)=SourceX(I,J)+diffusion

! ANNA start









! ANNA end

       IF(WindForce)THEN
          Tau_sx(I,J)=RHO_AW*Cdw*WindU2D(I,J) &
                *SQRT(WindU2D(I,J)*WindU2D(I,J) &
                      + WindV2D(I,J)*WindV2D(I,J))
          SourceX(I,J)=SourceX(I,J)+Jaco(I,J)*Tau_sx(I,J)

! ANNA start







! ANNA end

       ENDIF

       SourceY(I,J)=GRAV*(Eta(I,J))* &
          ((Depthy(I,J+1)*JL22y(I,J+1)-Depthy(I,J)*JL22y(I,J))/DY &
          +(Depthx(I+1,J)*JL12x(I+1,J)-Depthx(I,J)*JL12x(I,J))/DX)* &
          MASK(I,J) 


! ANNA start








! ANNA end

! friction
       Tau_by(I,J)=0.5_SP*Cd*(beta1*V(I,J) &
              +beta2*WaveUbottSC(I,J)*SIN(WaveAngleSC(I,J)*Deg2Rad))*Uu0

       SourceY(I,J)=SourceY(I,J)  &
          -Jaco(I,J)*Tau_by(I,J) &
                        ! Coriolis
                    -Coriolis(I,J)*Ubar(I,J) &
                       ! radiation stresses 
          +Jaco(I,J)*grav*WaveFySC(I,J)

! ANNA start














! ANNA end

! lateral mixing terms
       diffusion=-(nu_total(I,J)-Disp*D_22(I,J))  &
                    *NJL22*(2.0_SP*HL22*v_22+2.0_SP*v_2*HL22_2  &
                           +2.0_SP*HL12*v_12+2.0_SP*v_1*HL12_2 )  &
          -(nu_total(I,J)-Disp*D_22(I,J))  &
                    *NJL12*(2.0_SP*HL22*v_12+2.0_SP*v_2*HL22_1  &
                           +2.0_SP*HL12*v_11+2.0_SP*v_1*HL12_1 )  &
          -(nu_total(I,J)-Disp*D_11(I,J))  &
                    *NJL12*(HL12*v_22+v_2*HL12_2+HL22*v_12+v_1*HL11_2)  &
          -(nu_total(I,J)-Disp*D_22(I,J))  &
                    *NJL12*(HL22*u_22+u_2*HL22_2+HL12*u_12+u_1*HL12_2 )  &
          -(nu_total(I,J)-Disp*D_11(I,J))  &
                    *NJL11*(HL12*v_12+v_2*HL12_1+HL11*v_11+v_1*HL11_1)  &
          -(nu_total(I,J)-Disp*D_22(I,J))  &
                    *NJL11*(HL22*u_12+u_2*HL22_1+HL12*u_11+u_1*HL12_1 )
        SourceY(I,J)=SourceY(I,J)+diffusion  

! ANNA start









! ANNA end

       IF(WindForce)THEN
          Tau_sy(I,J)=RHO_AW*Cdw*WindV2D(I,J) &
                *SQRT(WindU2D(I,J)*WindU2D(I,J) &
                      + WindV2D(I,J)*WindV2D(I,J))
          SourceY(I,J)=SourceY(I,J)+Jaco(I,J)*Tau_sy(I,J)

! ANNA start







! ANNA end

       ENDIF

!     tmp4preview=SourceY

     ENDDO
     ENDDO

END SUBROUTINE SourceTerms

! --------------------------------------------------
!    This is subroutine to show statistics
!    called by
!        MAIN
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE STATISTICS
     USE GLOBAL
     IMPLICIT NONE

     REAL(SP)::MassVolume=ZERO,Energy=ZERO,MaxEta=ZERO,MinEta=ZERO, &
              MaxU=ZERO,MaxV=ZERO,Fr=ZERO,UTotal=ZERO,UTotalMax=ZERO



!
     MassVolume=ZERO
     Energy=ZERO
     UTotalMax=ZERO

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend

! Vol=SUM(Eta*dx*dy), reference is at z=0
! Energy=SUM(1/2*g*H^2*dx*dy+0.5*u^2*H*dx*dy)

       MassVolume=MassVolume+Eta(I,J)*DX*DY
       Energy=Energy+0.5_SP*H(I,J)*H(I,J)*GRAV*DX*DY &
             +0.5_SP*U(I,J)*U(I,J)*H(I,J)*DX*DY &
             +0.5_SP*V(I,J)*V(I,J)*H(I,J)*DX*DY      
     ENDDO
!      pause
     ENDDO
!     stop

     MaxEta=MAXVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MinEta=MINVAL(Eta(Ibeg:Iend,Jbeg:Jend))
     MaxU=MAXVAL(ABS(U(Ibeg:Iend,Jbeg:Jend)))
     MaxV=MAXVAL(ABS(V(Ibeg:Iend,Jbeg:Jend)))

! found Froude vs. max speed
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      IF(MASK(I,J)>ZERO)THEN
       Utotal=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
       IF(Utotal.gt.UtotalMax)THEN
         UtotalMax=Utotal
         Fr=SQRT(GRAV*Max(H(I,J),MinDepthfrc))
       ENDIF
      ENDIF
     ENDDO
     ENDDO
     IF(Fr==ZERO)Fr=SQRT(GRAV*MinDepthfrc)







! print screen
     WRITE(*,*)'----------------- STATISTICS ----------------'
     WRITE(*,*)' TIME        DT'
     WRITE(*,101) Time, DT
     WRITE(*,*)' MassVolume  Energy      MaxEta      MinEta      Max U       Max V '
     WRITE(*,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV
     WRITE(*,*)' MaxTotalU   PhaseS      Froude '
     WRITE(*,101) UTotalMax, Fr, UTotalMax/Fr
! print log file
     WRITE(3,*)'----------------- STATISTICS ----------------'
     WRITE(3,*)' TIME        DT'
     WRITE(3,101) Time, DT
     WRITE(3,*)' MassVolume  Energy      MaxEta      MinEta      Max U       Max V '
     WRITE(3,101) MassVolume,Energy,MaxEta,MinEta,MaxU,MaxV
     WRITE(3,*)' MaxTotalU   PhaseS      Froude '
     WRITE(3,101) UTotalMax, Fr, UTotalMax/Fr




101  FORMAT(6E12.4)

END SUBROUTINE STATISTICS

! ---------------------------------------------------
!    This is subroutine ESTIMATE_HUV 
!  for 3rd-order LK scheme 
!  called by
!      MAIN
!  call
!      GET_Eta_U_V_HU_HV
!    Last Update: 05/12/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE ESTIMATE_HUV(ISTEP)
     USE GLOBAL
     IMPLICIT NONE
     INTEGER,INTENT(IN)::ISTEP
     REAL(SP),PARAMETER::n_left=-1.0_SP,n_right=1.0_SP,n_bottom=-1.0_SP,n_top=1.0_SP
     REAL(SP)::F_left,F_right,F_bottom,F_top,WK_Source
     REAL(SP),DIMENSION(Ibeg:Iend,Jbeg:Jend)::R1,R2,R3
     INTEGER::kf

! MUSCL-Hancock, Zhou et al., p. 7

! solve eta
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend

      F_left=P(I,J)
      F_right=P(I+1,J)
      F_bottom=Q(I,J)
      F_top=Q(I,J+1)
         R1(I,J)=(-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                   -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom))&
                   /Max(SMALL, Jaco(I,J))








      Eta(I,J)=ALPHA(ISTEP)*Eta0(I,J)+BETA(ISTEP)*(Eta(I,J)+DT*R1(I,J))

     ENDDO
     ENDDO


! solve ubar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Fx(I,J)
      F_right=Fx(I+1,J)
      F_bottom=Fy(I,J)
      F_top=Fy(I,J+1)
      R2(I,J)=-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom) &
                        +SourceX(I,J)
      Ubar(I,J)=ALPHA(ISTEP)*Ubar0(I,J)+BETA(ISTEP)*(Ubar(I,J)+DT*R2(I,J))

! ANNA start











! ANNA end

     ENDDO
     ENDDO

! solve vbar
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      F_left=Gx(I,J)
      F_right=Gx(I+1,J)
      F_bottom=Gy(I,J)
      F_top=Gy(I,J+1)
      R3(I,J)=-1.0_SP/DX*(F_right*n_right+F_left*n_left) &
                       -1.0_SP/DY*(F_top*n_top+F_bottom*n_bottom) &
                       +SourceY(I,J)
      Vbar(I,J)=ALPHA(ISTEP)*Vbar0(I,J)+BETA(ISTEP)*(Vbar(I,J)+DT*R3(I,J))

! ANNA start


! ANNA end

     ENDDO
     ENDDO

     CALL GET_Eta_U_V_HU_HV

END SUBROUTINE ESTIMATE_HUV

! ---------------------------------------------------
!    This is subroutine to obtain Eta, u,v,hu,hv
!  called by
!      PREDICTOR
!      CORRECTOR
!      ESTIMATE_HUV (Lunge-Kutta)
!  use FroudeCap to Limit Froude<FroudeCap
!    Last Update: 12/06/2013 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE GET_Eta_U_V_HU_HV
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP)::Fr,Utotal,Utheta

! added 12/06/2013
      CALL UPDATE_MASK

! calculate etar, u and vetar, HU, HV
     H=Max(MinDepthFrc, Eta+Depth)

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend  

        U(I,J)=Ubar(I,J)/Max(H(I,J),MinDepthFrc)/Max(SMALL,Jaco(I,J))
        V(I,J)=Vbar(I,J)/Max(H(I,J),MinDepthFrc)/Max(SMALL,Jaco(I,J))




     ENDDO
     ENDDO   


     Uc = U*L11+V*L12
     Vc = U*L21+V*L22


     DO J=Jbeg,Jend
     DO I=Ibeg,Iend   
       IF(MASK(I,J)<1)THEN
        Ubar(I,J)=ZERO
        Vbar(I,J)=ZERO
        U(I,J)=ZERO
        V(I,J)=ZERO
        HU(I,J)=ZERO
        HV(I,J)=ZERO
       ELSE

        HU(I,J)=Max(H(I,J),MinDepthFrc)*Uc(I,J)
        HV(I,J)=Max(H(I,J),MinDepthFrc)*Vc(I,J)




! apply Froude cap
        Utotal=SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
        Fr=SQRT(GRAV*Max(H(I,J),MinDepthFrc))
        IF(Utotal/Fr.gt.FroudeCap)THEN
          Utheta=ATAN2(V(I,J),U(I,J))
          U(I,J)=FroudeCap*Fr*COS(Utheta)
          V(I,J)=FroudeCap*Fr*SIN(Utheta)

          HU(I,J)=(U(I,J)*L11(I,J)+V(I,J)*L12(I,J))*Max(H(I,J),MinDepthFrc)
          HV(I,J)=(U(I,J)*L21(I,J)+V(I,J)*L22(I,J))*Max(H(I,J),MinDepthFrc)




        ENDIF
! end Froude cap
       ENDIF
     ENDDO
     ENDDO

END SUBROUTINE GET_Eta_U_V_HU_HV

! ---------------------------------------------------
!    This is subroutine evaluate dt
!    Last Update: 12/09/2013 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE ESTIMATE_DT(M,N,DX,DY,U,V,H,MinDepthFrc,DT,CFL,TIME)
     USE PARAM
     USE PASS







     IMPLICIT NONE
     INTEGER,INTENT(IN)::M,N








     REAL(SP),DIMENSION(M,N),INTENT(IN)::DX,DY




     REAL(SP),INTENT(IN),DIMENSION(M,N)::U,V,H
     REAL(SP),INTENT(IN)::CFL,MinDepthFrc
     REAL(SP),INTENT(OUT)::DT
     REAL(SP),INTENT(INOUT)::TIME

     TMP3=LARGE
     DO J=1,N
     DO I=1,M
! x direction
      TMP1=ABS(U(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN


       TMP2=DX(I,J)/SMALL




      ELSE


       TMP2=DX(I,J)/TMP1




      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2
! y direction
      TMP1=ABS(V(I,J))+SQRT(GRAV*MAX(H(I,J),MinDepthFrc))
      IF(TMP1<SMALL)THEN


       TMP2=DY(I,J)/SMALL




      ELSE


       TMP2=DY(I,J)/TMP1




      ENDIF
      IF(TMP2<TMP3)TMP3=TMP2      
     ENDDO
     ENDDO







! added 12/09/2013







     DT=CFL*TMP3
! TEMP



     TIME=TIME+DT

     DT_SC=DT

END SUBROUTINE ESTIMATE_DT






! --------------------------------------------------
!    This is subroutine to damp waves using DHI type sponge layer 
!    variables
!    called by
!       MAIN
!    
!    Last Update: 10/27/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE SPONGE_DAMPING
     USE GLOBAL
     IMPLICIT NONE

     DO J=1,Nloc
     DO I=1,Mloc
      IF(MASK(I,J)>ZERO)THEN
       ETA(I,J)=ETA(I,J)/SPONGE(I,J)
      ENDIF
       U(I,J)=U(I,J)/SPONGE(I,J)
       V(I,J)=V(I,J)/SPONGE(I,J)
     ENDDO
     ENDDO

END SUBROUTINE SPONGE_DAMPING

! ------------------------------------------------
! This part is not subroutines
!  DEFINITIONS OF VARIABLES
! 
!    Last Update: 09/07/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
!
! Depth(): still water depth at element point
! DepthNode(): still water depth at node
! DepthX(): still water depth at x-interface
! DepthY(): still water depth at y-interface
! Eta():   surface elevation
! Eta0(): Eta at previous time level
!  for dry point, Eta() = MinDepth+Z()
! MASK(): 1 - wet
!         0 - dry
! MASK_STRUC(): 0 - permanent dry point
! MASK9: mask for itself and 8 elements around
! 
! U():  depth-averaged u or u at the reference level (u_alpha) at element
! V():  depth-averaged v or v at the reference level (v_alpha) at element
! HU(): (dep+eta)*u at element
! HV(): (dep+eta)*v at element
! P(): HU + dispersion at x-interface
! Q(): HV + dispersion at y-interface
! Fx(): F at x-interface
! Fy(): F at y-interface
! Gx(): G at x-interface
! Gy(): G at y-interface
! Ubar(:,:,:): Ubar
! Vbar(:,:,:): Vbar

! dispersion
! U1p(:,:): x-component of V1p
! V1p(:,:): y-component of V1p

! 
! EtaRxL(): Eta Left value at x-interface
! EtaRxR(): Eta Right value at x-interface
! EtaRyL(): Eta Left value at y-interface
! EtaRyR(): Eta Right value at y-interface
! HxL():   total depth  Left value at x-interface
! HxR():   total depth  Right value at x-interface
! HyL():   total depth  Left value at y-interface
! HyR():   total depth  Right value at y-interface

! HUxL(): HU Left value at x-interface
! HUxR(): HU Right value at x-interface
! HUyL(): HV Left value at y-interface
! HUyR(): HV Right value at y-interface

! PL(): HU + dispersion, Left value at x-interface
! PR(): HU + dispersion, Right value at x-interface
! QL(): HV + dispersion, Left value at y-interface
! QR(): HV + dispersion, Right value at y-interface

! FxL = HUxL*UxL + 1/2*g*(EtaRxL^2 + 2*EtaRxL*Depthx)
! FxR = HUxR*UxR + 1/2*g*(EtaRxR^2 + 2*EtaRxR*Depthx)
! FyL = HyL*UyL*VyL
! FyR = HyR*UyR*VyR

! GxL = HxL*UxL*VxL
! GxR = HxR*UxR*VxR
! GyL = HVyL*VyL + 1/2*g*(EtaRyL^2 + 2*EtaRyL*Depthy)
! GyR = HVyR*VyR + 1/2*g*(EtaRyR^2 + 2*EtaRyR*Depthy) 






