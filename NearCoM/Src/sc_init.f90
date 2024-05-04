! ---------------------------------------------------
!    This is subroutine to allocate variables
!  called by 
!        MAIN
!    Last Update: 05/07/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------

SUBROUTINE ALLOCATE_VARIABLES
     USE GLOBAL
     USE PASS




! allocate variables
     ALLOCATE(DelxU(Mloc,Nloc),DelxHU(Mloc,Nloc),DelxV(Mloc,Nloc),DelxEtar(Mloc,Nloc),&
              DelyU(Mloc,Nloc),DelyHV(Mloc,Nloc),DelyV(Mloc,Nloc),DelyEtar(Mloc,Nloc),&
              DelxHV(Mloc,Nloc),DelyHU(Mloc,Nloc), &
! U V HU H in x-direction
              UxL(Mloc1,Nloc),UxR(Mloc1,Nloc),VxL(Mloc1,Nloc),VxR(Mloc1,Nloc),&
              HUxL(Mloc1,Nloc),HUxR(Mloc1,Nloc),HVxL(Mloc1,Nloc),HVxR(Mloc1,Nloc), &
              HxL(Mloc1,Nloc),HxR(Mloc1,Nloc), &
! U V HV H in y-direction
              UyL(Mloc,Nloc1),UyR(Mloc,Nloc1),VyL(Mloc,Nloc1),VyR(Mloc,Nloc1),&
              HVyL(Mloc,Nloc1),HVyR(Mloc,Nloc1),HUyL(Mloc,Nloc1),HUyR(Mloc,Nloc1), &
              HyL(Mloc,Nloc1),HyR(Mloc,Nloc1), &
! P Q Eta, Fx, Fy
              PL(Mloc1,Nloc),PR(Mloc1,Nloc),QL(Mloc,Nloc1),QR(Mloc,Nloc1), &
              FxL(Mloc1,Nloc),FxR(Mloc1,Nloc),FyL(Mloc,Nloc1),FyR(Mloc,Nloc1), &
              GxL(Mloc1,Nloc),GxR(Mloc1,Nloc),GyL(Mloc,Nloc1),GyR(Mloc,Nloc1), &
              EtaRxL(Mloc1,Nloc),EtaRxR(Mloc1,Nloc), &
              EtaRyL(Mloc,Nloc1),EtaRyR(Mloc,Nloc1), &
! sponge
              SPONGE(Mloc,Nloc), &
! original variables at notes
              Fx(Mloc1,Nloc),Fy(Mloc,Nloc1),&
              U(Mloc,Nloc),V(Mloc,Nloc), HU(Mloc,Nloc),HV(Mloc,Nloc),&
              Gx(Mloc1,Nloc),Gy(Mloc,Nloc1), &
              P(Mloc1,Nloc),Q(Mloc,Nloc1), &
              SxL(Mloc1,Nloc),SxR(Mloc1,Nloc), &
              SyL(Mloc,Nloc1),SyR(Mloc,Nloc1),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc), &
              Coriolis(Mloc,Nloc) &
              )
! wave force
        IF(.NOT.ALLOCATED(WaveFxSC)) ALLOCATE(WaveFxSC(Mloc,Nloc))
        IF(.NOT.ALLOCATED(WaveFySC)) ALLOCATE(WaveFySC(Mloc,Nloc))

! others
! curvilinear

! fyshi changed X,Y dimension from Mglob to Mloc 04/02/2013
      ALLOCATE(X(Mloc,Nloc),Y(Mloc,Nloc), &
              Jaco(Mloc,Nloc),JL11x(Mloc1,Nloc),JL12x(Mloc1,Nloc), &
              JL21x(Mloc1,Nloc),JL22x(Mloc1,Nloc),JS1(Mloc,Nloc), &
              JS2(Mloc,Nloc),JacoX(Mloc1,Nloc),JacoY(Mloc,Nloc1), &
              JL11y(Mloc,Nloc1),JL12y(Mloc,Nloc1), &
              JL21y(Mloc,Nloc1),JL22y(Mloc,Nloc1), &
              Uc(Mloc,Nloc),Vc(Mloc,Nloc), &
              L11(Mloc,Nloc),L12(Mloc,Nloc),&
              L21(Mloc,Nloc),L22(Mloc,Nloc), &
              UnL(Mloc1,Nloc),UnR(Mloc1,Nloc), &
              VnL(Mloc,Nloc1),VnR(Mloc,Nloc1), &
              UbarXL(Mloc1,Nloc),UbarXR(Mloc1,Nloc),&
              UbarYL(Mloc,Nloc1),UbarYR(Mloc,Nloc1),&
              VbarXL(Mloc1,Nloc),VbarXR(Mloc1,Nloc),&
              VbarYL(Mloc,Nloc1),VbarYR(Mloc,Nloc1), &
              DelxUbar(Mloc,Nloc),DelxVbar(Mloc,Nloc),&
              DelyUbar(Mloc,Nloc),DelyVbar(Mloc,Nloc) &
              )

    
      ALLOCATE(D_11(Mloc,Nloc),D_22(Mloc,Nloc),D_12(Mloc,Nloc))

      ALLOCATE(Ubott(Mloc,Nloc),Vbott(Mloc,Nloc))

      ALLOCATE(Usurf(Mloc,Nloc),Vsurf(Mloc,Nloc))



      ALLOCATE(Depth(Mloc,Nloc),H(Mloc,Nloc),&
               Depthx(Mloc1,Nloc),Depthy(Mloc,Nloc1), &
               MASK(Mloc,Nloc), &
               MASK_STRUC(Mloc,Nloc),MASK9(Mloc,Nloc), &
               tmp4preview(Mloc,Nloc),Int2Flo(Mloc,Nloc)&
              )

      IF(FRC_MANNING_DATA)THEN
        ALLOCATE(CD_2D(Mloc,Nloc))
      ENDIF
! updating variables
      ALLOCATE(Eta(Mloc,Nloc),Eta0(Mloc,Nloc), &
               Ubar0(Mloc,Nloc),Vbar0(Mloc,Nloc),&
               Ubar(Mloc,Nloc),Vbar(Mloc,Nloc))

! mixing
      ALLOCATE(nu_total(Mloc,Nloc),nu_smg(Mloc,Nloc), &
               Tau_bx(Mloc,Nloc),Tau_by(Mloc,Nloc), &
               Tau_sx(Mloc,Nloc),Tau_sy(Mloc,Nloc))












END SUBROUTINE ALLOCATE_VARIABLES




! ----------------------------------------------------
!    This is subroutine initialization
!  called by 
!        MAIN
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIALIZATION
     USE GLOBAL



     USE Input_Util
     USE PASS
     IMPLICIT NONE
     CHARACTER(LEN=80) :: WHAT
     INTEGER :: RiverDirection,Ip,Jp
     REAL(SP),DIMENSION(Mglob,Nglob) :: DEPTH_noghost
     REAL(SP) :: TMP_L, TMP_R
     

     IF(.NOT.ALLOCATED(MASK_WC_INTERACT)) ALLOCATE(MASK_WC_INTERACT(Mglob,Nglob))
     MASK_WC_INTERACT = 1.0_SP
     DO J=1,Nglob
       DO I=1,WC_BOUND_WEST
         MASK_WC_INTERACT(I,J)=ZERO
       ENDDO
       DO I=Mglob-WC_BOUND_EAST+1,Mglob
         MASK_WC_INTERACT(I,J)=ZERO
       ENDDO
     ENDDO

     DO I=1,Mglob
       DO J=1,WC_BOUND_SOUTH
         MASK_WC_INTERACT(I,J)=ZERO
       ENDDO
       DO J=Nglob-WC_BOUND_NORTH+1,Nglob
         MASK_WC_INTERACT(I,J)=ZERO
       ENDDO
     ENDDO

! mask for bed update


! Coriolis parameter

     IF(CORI_CONSTANT)THEN

     Do J=1,Nloc
     Do I=1,Mloc
       Coriolis(I,J) = pi*SIN(LATITUDE*pi/180.0_SP) / 21600.0_SP   
     ENDDO
     ENDDO
   
     ELSE

! for serial
      OPEN(1,FILE=TRIM(LATITUDE_FILE))
       DO J=Jbeg,Jend
         READ(1,*)(Coriolis(I,J),I=Ibeg,Iend)
       ENDDO
      CLOSE(1)

! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        Coriolis(I,J)=Coriolis(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        Coriolis(I,J)=Coriolis(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        Coriolis(I,J)=Coriolis(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        Coriolis(I,J)=Coriolis(Iend,J)
       ENDDO
     ENDDO  

     Do J=1,Nloc
     Do I=1,Mloc
       Coriolis(I,J) = pi*SIN(Coriolis(I,J)*pi/180.0_SP) / 21600.0_SP   
     ENDDO
     ENDDO



     ENDIF















! set zeros






     Disp = 1.0_SP
     D_11 = ZERO
     D_12 = ZERO
     D_22 = ZERO
     nu_smg = ZERO
     DelxU=0.0_SP
     DelxHU=0.0_SP
     DelxV=0.0_SP
     DelxEtar=0.0_SP
     DelyU=0.0_SP
     DelyHV=0.0_SP
     DelyV=0.0_SP
     DelyEtar=0.0_SP
     DelxHV=0.0_SP
     DelyHU=0.0_SP
     UxL=0.0_SP
     UxR=0.0_SP
     VxL=0.0_SP
     VxR=0.0_SP
     HUxL=0.0_SP
     HUxR=0.0_SP
     HVxL=0.0_SP
     HVxR=0.0_SP
     HxL=0.0_SP
     HxR=0.0_SP
     UyL=0.0_SP
     UyR=0.0_SP
     VyL=0.0_SP
     VyR=0.0_SP
     HVyL=0.0_SP
     HVyR=0.0_SP
     HUyL=0.0_SP
     HUyR=0.0_SP
     HyL=0.0_SP
     HyR=0.0_SP

     PL=0.0_SP
     PR=0.0_SP
     QL=0.0_SP
     QR=0.0_SP
     FxL=0.0_SP
     FxR=0.0_SP
     FyL=0.0_SP
     FyR=0.0_SP
     GxL=0.0_SP
     GxR=0.0_SP
     GyL=0.0_SP
     GyR=0.0_SP
     SxL=0.0_SP
     SxR=0.0_SP
     SyL=0.0_SP
     SyR=0.0_SP
! original variables
     Ubar=0.0_SP
     Vbar=0.0_SP
     Ubar0=0.0_SP
     Vbar0=0.0_SP
     U=0.0_SP
     V=0.0_SP

     Uc=ZERO
     Vc=ZERO

     HU=0.0_SP
     HV=0.0_SP
     Fx=0.0_SP
     Fy=0.0_SP
     Gx=0.0_SP
     Gy=0.0_SP
     P=0.0_SP
     Q=0.0_SP
     Depth=10.0_SP
     H=0.0_SP
     Eta=0.0_SP
     SourceX=0.0_SP
     SourceY=0.0_SP
     PLOT_COUNT=0.0_SP
     PLOT_COUNT_STATION=0.0_SP
     HOTSTART_COUNT=ZERO
     MASK=1
     MASK_STRUC=1
     SCREEN_COUNT=ZERO
     SPONGE=1.0_SP


!  dx and dy become unit spacing if curvinear
     DX = 1.0_SP
     DY = 1.0_SP
! curvilinear stuff
     CALL JACOBIAN




     PLOT_COUNT=PLOT_INTV+SMALL
     PLOT_COUNT_STATION=PLOT_INTV_STATION+SMALL

     SCREEN_COUNT=SCREEN_INTV+SMALL
   
     DT=ZERO

! bathymetry

  IF(DEPTH_TYPE(1:3)=='DAT')THEN

! serial
     OPEN(1,FILE=TRIM(DEPTH_FILE))
       DO J=1,Nglob
        READ(1,*)(Depth_noghost(I,J),I=1,Mglob)
       ENDDO
     CLOSE(1)

     IF(PERIODIC_X)THEN
       DO J=1,NGlob
        TMP_L=Depth_noghost(MGlob-Num_Transit+1,J)
        TMP_R=Depth_noghost(Num_Transit,J)
        DO I=1,Num_Transit
          Depth_noghost(I,J)=TMP_L+(TMP_R-TMP_L)*(I+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO I=MGlob-Num_Transit,MGlob
          Depth_noghost(I,J)=TMP_L+(TMP_R-TMP_L)*(I-MGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

     IF(PERIODIC_Y)THEN
       DO I=1,MGlob
        TMP_L=Depth_noghost(I,NGlob-Num_Transit+1)
        TMP_R=Depth_noghost(I, Num_Transit)
        DO J=1,Num_Transit
          Depth_noghost(I,J)=TMP_L+(TMP_R-TMP_L)*(J+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO J=NGlob-Num_Transit,NGlob
          Depth_noghost(I,J)=TMP_L+(TMP_R-TMP_L)*(J-NGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

     DO J=Nghost+1,Nghost+Nglob
      DO I=Nghost+1,Nghost+Mglob
         Depth(I,J)=Depth_noghost(I-Nghost,J-Nghost)
      ENDDO
     ENDDO


  ENDIF

  IF(DEPTH_TYPE(1:3)=='FLA') THEN
    DO J=1,Nloc
     DO I=1,Mloc
      Depth(I,J) = Depth_FLat
     ENDDO
    ENDDO
  ENDIF



  IF(DEPTH_TYPE(1:3)=='SLO') THEN
    IF(INT(Xslp/DX)>Mloc.OR.INT(Xslp/DX)<1)THEN
      WRITE(*,*)'Please reset Xslp, STOP!'
      STOP
    ENDIF
    DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      Depth(I,J) = Depth_FLat
     ENDDO

     DO I=INT(Xslp/DX)+1+Ibeg,Iend
      Depth(I,J) = Depth_Flat-SLP*(I-(INT(Xslp/DX)+1+Ibeg))*DX
     ENDDO
    ENDDO
   ENDIF


! depth at ghost cells and re-construct depth at x,y-interfaces





! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        Depth(I,J)=Depth(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        Depth(I,J)=Depth(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        Depth(I,J)=Depth(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        Depth(I,J)=Depth(Iend,J)
       ENDDO
     ENDDO  



! ((((((((((((
  IF(FRC_MANNING_DATA)THEN





! serial
     OPEN(1,FILE=TRIM(FRC_FILE))
       DO J=1,Nglob
        READ(1,*)(CD_2D(I,J),I=1,Mglob)
       ENDDO
     CLOSE(1)

! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        CD_2D(I,J)=CD_2D(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        CD_2D(I,J)=CD_2D(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        CD_2D(I,J)=CD_2D(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        CD_2D(I,J)=CD_2D(Iend,J)
       ENDDO
     ENDDO  


  ENDIF

! ))))))))))))


! initial eta u and v for deforming
     IF(INI_UVZ) THEN
        CALL INITIAL_UVZ(Mloc,Nloc,Nghost,U_FILE,V_FILE,ETA_FILE,U,V,Eta)
     ENDIF 


! initial N wave
     IF(WaveMaker(1:6)=='N_WAVE') THEN
       CALL INITIAL_N_WAVE(Mloc,Nloc, DX,x1_Nwave,& 
          x2_Nwave,a0_Nwave,gamma_Nwave,dep_Nwave,U,V,Eta)
     ENDIF



! initial rectangular hump
     IF(WaveMaker(1:7)=='INI_REC') THEN
       CALL INITIAL_RECTANGULAR(Mloc,Nloc,Nghost,DX,DY,Xc,Yc,WID,AMP_SOLI, &
                      Eta)
     ENDIF



! initial wave surface
     IF(WaveMaker(1:7)=='INI_OTH') THEN
       CALL INITIAL_WAVE
     ENDIF



     IF(SPONGE_ON)THEN
       CALL CALCULATE_SPONGE(Mloc,Nloc,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)
     ENDIF

! water levels for stationary mode 



! get Eta and H
     H=MAX(Eta+Depth,MinDepthFrc)
     HU=H*U
     HV=H*V
  
     Ubar=HU
     Vbar=HV

     DO J=1,Nloc
     DO I=1,Mloc
      IF(Eta(I,J)<-DEPTH(I,J))THEN
       MASK(I,J)=0
       Eta(I,J)=MinDepth-Depth(I,J)
      ELSE
       MASK(I,J)=1
      ENDIF
     ENDDO
     ENDDO

! read obstacle structures
     IF(OBSTACLE)THEN
     IF(.not.check_exist(TRIM(OBSTACLE_FILE)))THEN
      !WRITE(*,*)TRIM(OBSTACLE_FILE), ' specified in input.txt but does not exist in folder.'
      STOP
     ENDIF
     OPEN(1,FILE=TRIM(OBSTACLE_FILE))
       DO J=Jbeg,Jend
        READ(1,*)(MASK_STRUC(I,J),I=Ibeg,Iend)
       ENDDO
     CLOSE(1)
     ENDIF

     DO J=1,Nloc
     DO I=1,Mloc
      IF(MASK_STRUC(I,J)==0)THEN
        Depth(I,J)=-LARGE
      ENDIF
     ENDDO
     ENDDO

! Tidal boundary conditions
     IF(ETA_CLAMPED)THEN
       OPEN(2,FILE=TRIM(TIDE_FILE))
         READ(2,*)WHAT

         READ(2,*) TideStartDate
         READ(2,*) NumConstituent
         ALLOCATE (TideFac(NumConstituent),Tideu0(NumConstituent))
         DO I = 1,NumConstituent
           READ(2,*) TideFac(I),Tideu0(I)
         ENDDO
         READ(2,*) NumEtaPoint
         ALLOCATE (I_eta_bd(NumEtaPoint),J_eta_bd(NumEtaPoint))
         ALLOCATE (TidePeriod(NumConstituent,NumEtaPoint), &
                   TideAmp(NumConstituent,NumEtaPoint),&
                   TidePha(NumConstituent,NumEtaPoint))
         DO I=1,NumEtaPoint
          READ(2,*)I_eta_bd(I),J_eta_bd(I)
          DO J=1,NumConstituent
            READ(2,*)TidePeriod(J,I),TideAmp(J,I),TidePha(J,I)
          ENDDO
         ENDDO
       CLOSE(2)
     ENDIF

! Tidal flux boundary conditions
     IF(FLUX_TIDE)THEN

       OPEN(2,FILE=TRIM(FLUX_TIDE_FILE))
         READ(2,*)WHAT

         READ(2,*) TideStartDate
         READ(2,*) NumConstituent
         READ(2,*) NumFluxPoint
         ALLOCATE (I_flux_bd(NumFluxPoint),J_flux_bd(NumFluxPoint), &
                   FLUX_ORIENTATION(NumFluxPoint), &
                   FluxAngle_Tide(NumFluxPoint) &
                   )
         ALLOCATE (FluxPeriod(NumConstituent,NumFluxPoint), &
                   FluxU_tide(NumConstituent,NumFluxPoint), &
                   FluxPha(NumConstituent,NumFluxPoint) &
                   )
         DO J=1,NumFluxPoint
          READ(2,*)I_flux_bd(J),J_flux_bd(J),FluxAngle_Tide(J),FLUX_ORIENTATION(J)

! !!!!!!!!!!
! adjust masks




         Ip=I_flux_bd(J)+Nghost
         Jp=J_flux_bd(J)+Nghost


         IF(FLux_Orientation(J)=='W')THEN
           IF (I_Flux_bd(J)/=1)THEN





             MASK_STRUC(Ip-1,Jp) = 0



           ENDIF 
         ELSEIF(Flux_Orientation(J)=='E')THEN
           IF (I_Flux_bd(J)/=Mglob)THEN





             MASK_STRUC(Ip+1,Jp) = 0



           ENDIF 
         ELSEIF(Flux_Orientation(J)=='S')THEN
           IF (J_Flux_bd(J)/=1)THEN





             MASK_STRUC(Ip,Jp-1) = 0



           ENDIF 
         ELSEIF(Flux_Orientation(J)=='N')THEN
           IF (J_Flux_bd(J)/=Nglob)THEN





             MASK_STRUC(Ip,Jp+1) = 0



           ENDIF 
         ENDIF

! adjust over 

! !!!!!!!!!
          DO I=1,NumConstituent
            READ(2,*)FluxPeriod(I,J),FluxU_tide(I,J),FluxPha(I,J)



          ENDDO
         ENDDO
       CLOSE(2)
     ENDIF

! wind
! wind force ! now constant wind, can use WindU2D for 2D applications
      IF(WindForce)THEN
        ALLOCATE(WindU2D(Mloc,Nloc),WindV2D(Mloc,Nloc))
       OPEN(2,FILE=TRIM(WIND_FILE))
         READ(2,*)WHAT
         READ(2,*)NumTimeWindData
         ALLOCATE (TimeWind(NumTimeWindData),WU(NumTimeWindData),&
                WV(NumTimeWindData))
         DO I=1,NumTimeWindData
           READ(2,*,END=111)TimeWind(I),WU(I),WV(I)
         ENDDO
111    CONTINUE
       CLOSE(2)
      ENDIF


! river flux boundary conditions
     IF(FLUX_CLAMPED)THEN
       OPEN(2,FILE=TRIM(FLUX_FILE))
         READ(2,*)WHAT
         READ(2,*)NumTimeData,NumFluxPoint
         ALLOCATE (FluxU_bd(NumFluxPoint,NumTimeData), &
                    FluxAngle(NumFluxPoint,NumTimeData))
         ALLOCATE (I_flux_bd(NumFluxPoint),&
                   J_flux_bd(NumFluxPoint),&
                   RIVER_ORIENTATION(NumFluxPoint) )
         ALLOCATE (TimeFlux(NumTimeData))

         DO J=1,NumFluxPoint
           READ(2,*)I_flux_bd(J),J_flux_bd(J),RIVER_ORIENTATION(J)

! adjust masks




         Ip=I_flux_bd(J)+Nghost
         Jp=J_flux_bd(J)+Nghost


         IF(River_Orientation(J)=='W')THEN
           IF (I_Flux_bd(J)/=1)THEN





             MASK_STRUC(Ip-1,Jp) = 0



           ENDIF 
         ELSEIF(River_Orientation(J)=='E')THEN
           IF (I_Flux_bd(J)/=Mglob)THEN





             MASK_STRUC(Ip+1,Jp) = 0



           ENDIF 
         ELSEIF(River_Orientation(J)=='S')THEN
           IF (J_Flux_bd(J)/=1)THEN





             MASK_STRUC(Ip,Jp-1) = 0



           ENDIF 
         ELSEIF(River_Orientation(J)=='N')THEN
           IF (J_Flux_bd(J)/=Nglob)THEN





             MASK_STRUC(Ip,Jp+1) = 0



           ENDIF 
         ENDIF

! adjust over 
           DO I=1,NumTimeData
             READ(2,*) TimeFlux(I), FluxU_bd(J,I),FluxAngle(J,I)
           ENDDO

         ENDDO ! end num flux point

      ENDIF ! end flux_clamped

! masks 

     MASK=MASK*MASK_STRUC

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
      MASK9(I,J)=MASK(I,J)*MASK(I-1,J)*MASK(I+1,J)  &
                *MASK(I+1,J+1)*MASK(I,J+1)*MASK(I-1,J+1) &
                *MASK(I+1,J-1)*MASK(I,J-1)*MASK(I-1,J-1) 
     ENDDO
     ENDDO

! re-construct Depth

     DO J=1,Nloc
     DO I=2,Mloc
      DepthX(I,J)=0.5_SP*(Depth(I-1,J)+Depth(I,J))
     ENDDO
     ENDDO
     DO J=1,Nloc
      DepthX(1,J)=0.5_SP*(3.0_SP*Depth(1,J)-Depth(2,J))
      DepthX(Mloc1,J)=0.5_SP*(3.0_SP*Depth(Mloc,J)-Depth(Mloc-1,J))
     ENDDO

     DO J=2,Nloc
     DO I=1,Mloc
      DepthY(I,J)=0.5_SP*(Depth(I,J-1)+Depth(I,J))
     ENDDO
     ENDDO
     DO I=1,Mloc
      DepthY(I,1)=0.5_SP*(3.0_SP*Depth(I,1)-Depth(I,2))
      DepthY(I,Nloc1)=0.5_SP*(3.0_SP*Depth(I,Nloc)-Depth(I,Nloc-1))
     ENDDO
      
!   don't need to re-construct depth in terms of using well-balanced scheme
!     DO J=1,Nloc
!     DO I=1,Mloc
!       Depth(I,J)=0.25_SP*(Depthx(I,J)+Depthx(I+1,J)+Depthy(I,J)+Depthy(I,J+1))
!     ENDDO
!     ENDDO

! deal with masks

      DO J=2,Nloc-1
      DO I=2,Mloc-1
        IF(MASK(I,J)<1)THEN
         DepthX(I,J)=Depth(I-1,J)
         DepthX(I+1,J)=Depth(I+1,J)
         DepthY(I,J)=Depth(I,J-1)
         DepthY(I,J+1)=Depth(I,J+1)
        ENDIF
      ENDDO
      ENDDO




 
       OPEN(11,FILE=TRIM(COUPLING_FILE))
         READ(11,*)  ! title
         READ(11,*)  ! boundary info
! boundary basic info including point number of coupling, start point, etc
! east
         READ(11,*)  ! east
         READ(11,*) N_COUPLING_EAST,J_START_EAST
! west 
         READ(11,*)  ! west
         READ(11,*) N_COUPLING_WEST,J_START_WEST
! south 
         READ(11,*)  ! south
         READ(11,*) N_COUPLING_SOUTH,I_START_SOUTH
! north 
         READ(11,*)  ! north
         READ(11,*) N_COUPLING_NORTH,I_START_NORTH

! read time and variable at the first level

         READ(11,*) ! time start title
         READ(11,*) TIME_COUPLING_1 
! initialize time_2
         TIME_COUPLING_2 = TIME_COUPLING_1

! east
         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,2),&
               V_COUPLING_EAST(N_COUPLING_EAST,2),&
               Z_COUPLING_EAST(N_COUPLING_EAST,2))
             READ(11,*)   ! east
             READ(11,119)(U_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)(V_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
!   initialize first step
             U_COUPLING_EAST(:,1)=U_COUPLING_EAST(:,2)
             V_COUPLING_EAST(:,1)=V_COUPLING_EAST(:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_east
119      FORMAT(5E16.6)

! west
         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,2),&
               V_COUPLING_WEST(N_COUPLING_WEST,2),&
               Z_COUPLING_WEST(N_COUPLING_WEST,2))
             READ(11,*)   ! west
             READ(11,119)(U_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)(V_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
!   initialize first step
             U_COUPLING_WEST(:,1)=U_COUPLING_WEST(:,2)
             V_COUPLING_WEST(:,1)=V_COUPLING_WEST(:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_west
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,2))
             READ(11,*)   ! south
             READ(11,119)(U_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)(V_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
!   initialize first step
             U_COUPLING_SOUTH(:,1)=U_COUPLING_SOUTH(:,2)
             V_COUPLING_SOUTH(:,1)=V_COUPLING_SOUTH(:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_south
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,2),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,2),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,2))
             READ(11,*)   ! north
             READ(11,119)(U_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)(V_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)

!   initialize first step
             U_COUPLING_NORTH(:,1)=U_COUPLING_NORTH(:,2)
             V_COUPLING_NORTH(:,1)=V_COUPLING_NORTH(:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_north


! specify boundary start points

! west boundary
   IF(N_COUPLING_WEST>0)THEN

      Kstart_WEST=J_START_WEST+Nghost
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1
      Kshift_WEST = -(Kstart_WEST-Nghost)+1
      IN_DOMAIN_WEST = .TRUE.


   ENDIF

! east boundary
   IF(N_COUPLING_EAST>0)THEN

      Kstart_EAST=J_START_EAST+Nghost
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1
      Kshift_EAST = -(Kstart_EAST-Nghost)+1
      IN_DOMAIN_EAST = .TRUE.

    ENDIF

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN

      Kstart_SOUTH=I_START_SOUTH+Nghost
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1
      Kshift_SOUTH = -(Kstart_SOUTH-Nghost)+1
      IN_DOMAIN_SOUTH = .TRUE.

   ENDIF

! north boundary
   IF(N_COUPLING_NORTH>0)THEN

      Kstart_NORTH=I_START_NORTH+Nghost
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1
      Kshift_NORTH = -(Kstart_NORTH-Nghost)+1
      IN_DOMAIN_NORTH = .TRUE.

   ENDIF


! end coupling


END SUBROUTINE INITIALIZATION

! ----------------------------------------------------
!    This is subroutine to index for MPI
!  called by 
!        MAIN
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------

SUBROUTINE INDEX
    USE GLOBAL
    IMPLICIT NONE

! TEMP

  px=1
  py=1


! check if Mloc and Nloc are the same for every processor









    Mloc=Mglob/px+2*Nghost
    Nloc=Nglob/py+2*Nghost
    Mloc1=Mloc+1
    Nloc1=Nloc+1

    Ibeg=Nghost+1
    Iend=Mloc-Nghost
    Iend1=Mloc1-Nghost
    Jbeg=Nghost+1
    Jend=Nloc-Nghost
    Jend1=Nloc1-Nghost

!    CALL SWEXITMPI
!    STOP

END SUBROUTINE INDEX

! ---------------------------------------------------
!    This is subroutine of sponge layer to get SPONGE
!  called by
!      INITIALIZATION
!    Last Update: 10/27/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE CALCULATE_SPONGE(M,N,Nghost,DX,DY,&
                            Sponge_west_width,Sponge_east_width,&
                            Sponge_south_width,Sponge_north_width, &
                            R_sponge,A_sponge,SPONGE)
     USE PARAM



     IMPLICIT NONE
     INTEGER, INTENT(IN)::M,N,Nghost
     REAL(SP),INTENT(IN)::DX,DY

     REAL(SP),INTENT(IN) :: &
                          Sponge_west_width,Sponge_east_width,&
                          Sponge_south_width,Sponge_north_width, &
                          R_sponge,A_sponge
     REAL(SP),DIMENSION(M,N),INTENT(INOUT)::SPONGE
     REAL(SP)::ri,lim
     INTEGER::Iwidth

! west



     IF(Sponge_west_width>ZERO)THEN

     Iwidth=INT(Sponge_west_width/DX)+Nghost

     DO J=1,N
     DO I=1,Iwidth
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
       ri=R_sponge**(50*(I-1)/(Iwidth-1))
       SPONGE(I,J)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF



! east



     IF(Sponge_east_width>ZERO)THEN

     Iwidth=INT(Sponge_east_width/DX)+Nghost

     DO J=1,N
     DO I=M-Iwidth+1,M
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
       ri=R_sponge**(50*(M-I)/(Iwidth-1))
       SPONGE(I,J)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF




! south



     IF(Sponge_south_width>ZERO)THEN

     Iwidth=INT(Sponge_south_width/DY)+Nghost

     DO I=1,M
     DO J=1,Iwidth
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
       ri=R_sponge**(50*(J-1)/(Iwidth-1))
       SPONGE(I,J)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF




! north



     IF(Sponge_north_width>ZERO)THEN

     Iwidth=INT(Sponge_north_width/DY)+Nghost

     DO I=1,M
     DO J=N-Iwidth+1,N
       IF(SPONGE(I,J)>1.0_SP)THEN
         lim=SPONGE(I,J)
       ELSE
         lim=1.0_SP
       ENDIF
       ri=R_sponge**(50*(N-J)/(Iwidth-1))
       SPONGE(I,J)=MAX(A_Sponge**ri,lim)
     ENDDO
     ENDDO
     ENDIF




END SUBROUTINE CALCULATE_SPONGE


! ---------------------------------------------------
!    This is subroutine to calculate Jacobian etc
!  called by
!      INITIALIZATION
!    Last Update: 05/25/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE JACOBIAN
      USE GLOBAL
      IMPLICIT NONE





      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Xglob,Yglob, &
              Jglob,JL11gx,JL12gx,JL21gx,JL22gx, &
              JL11gy,JL12gy,JL21gy,JL22gy,&
              JS1g,JS2g,JglobX, L11g,L12g,L21g,L22g,&
              JglobY

      ALLOCATE (Xglob(MGlob+2*Nghost,NGlob+2*Nghost), &
                Yglob(MGlob+2*Nghost,NGlob+2*Nghost), &
                Jglob(MGlob+2*Nghost,NGlob+2*Nghost), &
                JglobX(MGlob+2*Nghost+1,NGlob+2*Nghost), &
                JglobY(MGlob+2*Nghost,NGlob+2*Nghost+1), &
                JL11gx(MGlob+2*Nghost+1,NGlob+2*Nghost), &
                JL12gx(MGlob+2*Nghost+1,NGlob+2*Nghost), &
                JL21gx(MGlob+2*Nghost+1,NGlob+2*Nghost), &
                JL22gx(MGlob+2*Nghost+1,NGlob+2*Nghost), &
                JL11gy(MGlob+2*Nghost,NGlob+2*Nghost+1), &
                JL12gy(MGlob+2*Nghost,NGlob+2*Nghost+1), &
                JL21gy(MGlob+2*Nghost,NGlob+2*Nghost+1), &
                JL22gy(MGlob+2*Nghost,NGlob+2*Nghost+1), &
                JS1g(MGlob+2*Nghost,NGlob+2*Nghost), &
                JS2g(MGlob+2*Nghost,NGlob+2*Nghost),&
                L11g(MGlob+2*Nghost,NGlob+2*Nghost), &
                L12g(MGlob+2*Nghost,NGlob+2*Nghost), &
                L21g(MGlob+2*Nghost,NGlob+2*Nghost), &
                L22g(MGlob+2*Nghost,NGlob+2*Nghost) &
                  )

!# if defined (PARALLEL) every processor may do this
!     if (myid.eq.0) then
!# endif
        OPEN(1,FILE=TRIM(X_FILE))
        DO J=Nghost+1,NGlob+NGhost
           READ(1,*)(Xglob(I,J),I=Nghost+1,MGlob+Nghost)
        ENDDO
        CLOSE(1)

        OPEN(1,FILE=TRIM(Y_FILE))
        DO J=Nghost+1,NGlob+NGhost
           READ(1,*)(Yglob(I,J),I=Nghost+1,MGlob+Nghost)
        ENDDO
        CLOSE(1)

! J=X_xiY_eta - X_eta Y_xi
! JL11= Y_eta,  JL12 = - X_eta
! JL21 = - Y_xi,JL22= X_xi 
! where Lij=L^i_j
 
      DO J=Nghost+2,Nglob+Nghost-1
      DO I=Nghost+2,Mglob+Nghost-1
        Jglob(I,J) = 0.5_SP*(Xglob(I+1,J)-Xglob(I-1,J)) &
                    *0.5_SP*(Yglob(I,J+1)-Yglob(I,J-1)) &
                   -0.5_SP*(Xglob(I,J+1)-Xglob(I,J-1)) &
                    *0.5_SP*(Yglob(I+1,J)-Yglob(I-1,J))
        JS1g(I,J) = SQRT(0.25_SP*(Xglob(I+1,J)-Xglob(I-1,J))**2  &
                       +0.25_SP*(Yglob(I+1,J)-Yglob(I-1,J))**2)
        JS2g(I,J) = SQRT(0.25_SP*(Xglob(I,J+1)-Xglob(I,J-1))**2  &
                       +0.25_SP*(Yglob(I,J+1)-Yglob(I,J-1))**2)
      ENDDO
      ENDDO

! x point
      DO J=Nghost+2,Nglob+Nghost-1
      DO I=Nghost+2,Mglob+Nghost
        tmp1=0.5_SP*(Yglob(I,J+1)+Yglob(I-1,J+1))
        tmp2=0.5_SP*(Yglob(I,J-1)+Yglob(I-1,J-1))
        tmp3=0.5_SP*(Xglob(I,J+1)+Xglob(I-1,J+1))
        tmp4=0.5_SP*(Xglob(I,J-1)+Xglob(I-1,J-1))
        JglobX(I,J)=1.0_SP*(Xglob(I,J)-Xglob(I-1,J)) &
                    *0.5_SP*(tmp1-tmp2) &
                   -0.5_SP*(tmp3-tmp4) &
                    *1.0_SP*(Yglob(I,J)-Yglob(I-1,J))
        JL11gx(I,J) = 0.5_SP*(tmp1-tmp2)
        JL12gx(I,J) =-0.5_SP*(tmp3-tmp4)
        JL21gx(I,J) =-1.0_SP*(Yglob(I,J)-Yglob(I-1,J))
        JL22gx(I,J) = 1.0_SP*(Xglob(I,J)-Xglob(I-1,J))
      ENDDO
      ENDDO

! y point
      DO J=Nghost+2,Nglob+Nghost
      DO I=Nghost+2,Mglob+Nghost-1
        tmp1=0.5_SP*(Xglob(I+1,J)+Xglob(I+1,J-1))
        tmp2=0.5_SP*(Xglob(I-1,J)+Xglob(I-1,J-1))
        tmp3=0.5_SP*(Yglob(I+1,J)+Yglob(I+1,J-1))
        tmp4=0.5_SP*(Yglob(I-1,J)+Yglob(I-1,J-1))
        JglobY(I,J)=1.0_SP*(Yglob(I,J)-Yglob(I,J-1)) &
                    *0.5_SP*(tmp1-tmp2) &
                   -0.5_SP*(tmp3-tmp4) &
                    *1.0_SP*(Xglob(I,J)-Xglob(I,J-1))
        JL22gy(I,J) = 0.5_SP*(tmp1-tmp2)
        JL21gy(I,J) =-0.5_SP*(tmp3-tmp4)
        JL11gy(I,J)= 1.0_SP*(Yglob(I,J)-Yglob(I,J-1))
        JL12gy(I,J)=-1.0_SP*(Xglob(I,J)-Xglob(I,J-1))
      ENDDO
      ENDDO

! ghost cells 
! centroid
        DO I=Nghost+2,MGlob+Nghost-1
           DO J=1,Nghost+1
              Jglob(I,J)=Jglob(I,Nghost+2)
              JS1g(I,J)=JS1g(I,Nghost+2)
              JS2g(I,J)=JS2g(I,Nghost+2)
           ENDDO
           DO J=NGlob+Nghost,NGlob+2*Nghost
              Jglob(I,J)=Jglob(I,NGlob+Nghost-1)
              JS1g(I,J)=JS1g(I,NGlob+Nghost-1)
              JS2g(I,J)=JS2g(I,NGlob+Nghost-1)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost
           DO I=1,Nghost+1
              Jglob(I,J)=Jglob(Nghost+2,J)
              JS1g(I,J)=JS1g(Nghost+2,J)
              JS2g(I,J)=JS2g(Nghost+2,J)
           ENDDO
           DO I=MGlob+Nghost,MGlob+2*Nghost
              Jglob(I,J)=Jglob(MGlob+Nghost-1,J)
              JS1g(I,J)=JS1g(MGlob+Nghost-1,J)
              JS2g(I,J)=JS2g(MGlob+Nghost-1,J)
           ENDDO
        ENDDO
! x point
        DO I=Nghost+2,MGlob+Nghost
           DO J=1,Nghost+1
              JL11gx(I,J)=JL11gx(I,Nghost+2)
              JL12gx(I,J)=JL12gx(I,Nghost+2)
              JL21gx(I,J)=JL21gx(I,Nghost+2)
              JL22gx(I,J)=JL22gx(I,Nghost+2)
              JglobX(I,J)=JglobX(I,Nghost+2)
           ENDDO
           DO J=NGlob+Nghost,NGlob+2*Nghost
              JL11gx(I,J)=JL11gx(I,NGlob+Nghost-1)
              JL12gx(I,J)=JL12gx(I,NGlob+Nghost-1)
              JL21gx(I,J)=JL21gx(I,NGlob+Nghost-1)
              JL22gx(I,J)=JL22gx(I,NGlob+Nghost-1)
              JglobX(I,J)=JglobX(I,NGlob+Nghost-1)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost
           DO I=1,Nghost+1
              JL11gx(I,J)=JL11gx(Nghost+2,J)
              JL12gx(I,J)=JL12gx(Nghost+2,J)
              JL21gx(I,J)=JL21gx(Nghost+2,J)
              JL22gx(I,J)=JL22gx(Nghost+2,J)
              JglobX(I,J)=JglobX(Nghost+2,J)
           ENDDO
           DO I=MGlob+Nghost+1,MGlob+2*Nghost+1
              JL11gx(I,J)=JL11gx(MGlob+Nghost,J)
              JL12gx(I,J)=JL12gx(MGlob+Nghost,J)
              JL21gx(I,J)=JL21gx(MGlob+Nghost,J)
              JL22gx(I,J)=JL22gx(MGlob+Nghost,J)
              JglobX(I,J)=JglobX(MGlob+Nghost,J)
           ENDDO
        ENDDO
! y point
        DO I=Nghost+2,MGlob+Nghost-1
           DO J=1,Nghost+1
              JL21gy(I,J)=JL21gy(I,Nghost+2)
              JL22gy(I,J)=JL22gy(I,Nghost+2)
              JL11gy(I,J)=JL11gy(I,Nghost+2)
              JL12gy(I,J)=JL12gy(I,Nghost+2)
              JglobY(I,J)=JglobY(I,Nghost+2)
           ENDDO
           DO J=NGlob+Nghost+1,NGlob+2*Nghost+1
              JL21gy(I,J)=JL21gy(I,NGlob+Nghost)
              JL22gy(I,J)=JL22gy(I,NGlob+Nghost)
              JL11gy(I,J)=JL11gy(I,NGlob+Nghost)
              JL12gy(I,J)=JL12gy(I,NGlob+Nghost)
              JglobY(I,J)=JglobY(I,NGlob+Nghost)
           ENDDO
        ENDDO
        DO J=1,NGlob+2*Nghost+1
           DO I=1,Nghost+1
              JL21gy(I,J)=JL21gy(Nghost+2,J)
              JL22gy(I,J)=JL22gy(Nghost+2,J)
              JL11gy(I,J)=JL11gy(Nghost+2,J)
              JL12gy(I,J)=JL12gy(Nghost+2,J)
              JglobY(I,J)=JglobY(Nghost+2,J)
           ENDDO
           DO I=MGlob+Nghost,MGlob+2*Nghost
              JL21gy(I,J)=JL21gy(MGlob+Nghost-1,J)
              JL22gy(I,J)=JL22gy(MGlob+Nghost-1,J)
              JL11gy(I,J)=JL11gy(MGlob+Nghost-1,J)
              JL12gy(I,J)=JL12gy(MGlob+Nghost-1,J)
              JglobY(I,J)=JglobY(MGlob+Nghost-1,J)
           ENDDO
        ENDDO

        DO J=1,NGlob+2*Nghost
        DO I=1,MGlob+2*Nghost
         L11g(I,J)=(JL11gx(I,J)+JL11gx(I+1,J)+JL11gy(I,J)+JL11gy(I,J+1)) &
                     *0.25_SP/MAX(SMALL,Jglob(I,J))
         L12g(I,J)=(JL12gx(I,J)+JL12gx(I+1,J)+JL12gy(I,J)+JL12gy(I,J+1)) &
                     *0.25_SP/MAX(SMALL,Jglob(I,J))
         L21g(I,J)=(JL21gx(I,J)+JL21gx(I+1,J)+JL21gy(I,J)+JL21gy(I,J+1)) &
                     *0.25_SP/MAX(SMALL,Jglob(I,J))
         L22g(I,J)=(JL22gx(I,J)+JL22gx(I+1,J)+JL22gy(I,J)+JL22gy(I,J+1)) &
                     *0.25_SP/MAX(SMALL,Jglob(I,J))
        ENDDO
        ENDDO
     
!# if defined (PARALLEL)  every processor may do this
!     endif ! end myid=0
!# endif


      X = Xglob
      Y = Yglob

      Jaco = Jglob
      JL11x = JL11gx
      JL12x = JL12gx
      JL21x = JL21gx
      JL22x = JL22gx
      JL11y = JL11gy
      JL12y = JL12gy
      JL21y = JL21gy
      JL22y = JL22gy
      L11 = L11g
      L12 = L12g
      L21 = L21g
      L22 = L22g
      JS1  = JS1g
      JS2  = JS2g
      JacoX = JglobX
      JacoY = JglobY


      DEALLOCATE (Xglob, &
                Yglob, &
                Jglob, &
                JglobX, &
                JglobY, &
                JL11gx, &
                JL12gx, &
                JL21gx, &
                JL22gx, &
                JL11gy, &
                JL12gy, &
                JL21gy, &
                JL22gy, &
                L11g, &
                L12g, &
                L21g, &
                L22g, &
                JS1g,JS2g )

END SUBROUTINE JACOBIAN



! ---------------------------------------------------
!    This is subroutine of given initial u v and eta 
!  called by
!      INITIALIZATION
!    Last Update: 02/03/2011 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIAL_UVZ(M,N,Nghost,U_FILE,V_FILE,ETA_FILE,U,V,ETA)
      USE PARAM
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,Nghost
      REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U,V,ETA
      CHARACTER(LEN=80),INTENT(IN) :: ETA_FILE
      CHARACTER(LEN=80),INTENT(IN) :: U_FILE
      CHARACTER(LEN=80),INTENT(IN) :: V_FILE







      OPEN(1,FILE=TRIM(U_FILE))
       DO J=Nghost+1,N-Nghost
        READ(1,*)(U(I,J),I=Nghost+1,M-Nghost)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(V_FILE))
       DO J=Nghost+1,N-Nghost
        READ(1,*)(V(I,J),I=Nghost+1,M-Nghost)
       ENDDO
      CLOSE(1)

      OPEN(1,FILE=TRIM(ETA_FILE))
       DO J=Nghost+1,N-Nghost
        READ(1,*)(ETA(I,J),I=Nghost+1,M-Nghost)
       ENDDO
      CLOSE(1)


END SUBROUTINE INITIAL_UVZ


! --------------------------------------------------
!    This is subroutine to provide initial wave condition
!    it can be specified in input.txt giving 'INI'
!    called by
!       MAIN
!    Last Update: 10/01/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIAL_WAVE
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP),Dimension(Mloc,Nloc)::XX,YY
     REAL(SP):: sigma,x_c,y_c,a,r
     INTEGER :: ii1,ii2

! the initial wave domain includes ghost cells

     do j=1,Nloc
     do i=1,Mloc
       xx(i,j)=(i-1.)*dx
       yy(i,j)=(j-1.)*dy
     enddo
     enddo

     sigma=0.5
     x_c=10.0
     y_c=10.0	
     a=0.1

     do j=1,Nloc
     do i=1,Mloc
     goto 200
! box
     tmp1=21+Nghost
     tmp2=31+Nghost
     if(i>tmp1.and.i<tmp2.and.j>tmp1.and.j<tmp2)then
      Eta(i,j)=1.0_SP
     else
      Eta(i,j)=zero
     endif

200  continue

     goto 100
! dam break
        if(i<100+Nghost)then
         Eta(i,j)=5.0_SP
        else
         Eta(i,j)=0.0_SP
        endif
100   continue
! gausian
!         r=sqrt((xx(i,j)-x_c)**2+(yy(i,j)-y_c)**2)
!         Eta(i,j)=a*exp(-r**2/2./sigma**2)
     enddo
     enddo

! alongshore crest
     goto 213
     ii1=21+Nghost
     ii2=25+Nghost
     DO J=1,Nloc
       DO I=ii1,ii2
         Eta(I,J)=1.0_SP
       ENDDO
     ENDDO
213  continue

END SUBROUTINE INITIAL_WAVE



! --------------------------------------------------
!    This is subroutine to provide initial rectangular hump 
!    it can be specified in input.txt giving 'INI_REC'
!    called by
!       - INITIALIZATION
!    Last Update: 10/11/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIAL_RECTANGULAR(M,N,Nghost,DX,DY,Xc,Yc,WID,AMP,Eta)
     USE PARAM



     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N,Nghost
     REAL(SP),INTENT(IN) :: DX,DY,Xc,Yc,WID,AMP
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: Eta
     INTEGER :: Il,Jl,Ir,Jr

     Eta = ZERO


     Il=Xc/DX+1 - WID/DX + Nghost
     Ir=Xc/DX+1 + WID/DX + Nghost
     Jl=Yc/DY+1 - WID/DY + Nghost
     Jr=Yc/DY+1 + WID/DY + Nghost

     Eta(Il:Ir,Jl:Jr) = AMP

 
END SUBROUTINE INITIAL_RECTANGULAR



! --------------------------------------------------
!    This is subroutine to provide initial N wave solution 
!    it can be specified in input.txt giving 'N_WAVE'
!    called by
!       - INITIAL_N_WAVE
!    Last Update: 11/30/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIAL_N_WAVE(M,N,DX,x1,x2,a0,gamma,dep,U,V,Eta)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N
     REAL(SP),INTENT(IN) :: x1,x2,a0,gamma,dep
     REAL(SP),INTENT(IN) :: DX

     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U,V,Eta
     IF(dep==ZERO)THEN
      !WRITE(*,*)'depth for Nwave should not be 0, stop!'
      STOP
     ENDIF
     DO J=1,N
     DO I=1,M
       tmp1=(I-1)*DX
       Eta(I,J)=a0*(tmp1-x2)/COSH(gamma*(tmp1-x1))
       U(I,J) = SQRT(GRAV/dep)*Eta(I,J)
       V(I,J) = ZERO
     ENDDO
     ENDDO

END SUBROUTINE INITIAL_N_WAVE



! --------------------------------------------------
!    This is subroutine to provide initial solitary wave solution 
!    it can be specified in input.txt giving 'SOL'
!    called by
!       - INITIAL_WAVE
!    call
!       - SOLITARY_SOLUTION
!    Last Update: 10/01/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE INITIAL_SOLITARY_WAVE(M,N,DX,Xwavemaker,&
           AMP_Soli,Dep_Soli,Beta,U,V,Eta)
     USE PARAM
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: M,N
     REAL(SP),INTENT(IN) :: DX 
     REAL(SP),INTENT(IN) :: Dep_Soli,Xwavemaker,Beta
     REAL(SP),DIMENSION(M,N),INTENT(OUT) :: U,V,Eta
     REAL(SP),INTENT(IN) :: AMP_Soli
     REAL(SP) :: A,B,A1,A2,alpha,SC,Cph,C_ph

!     alpha=0.5_SP*Beta*Beta+Beta
     alpha=-0.39

      CALL SUB_SLTRY(AMP_Soli,Dep_Soli,alpha,Cph,B,A1,A2,A,C_ph)

     DO J=1,N
     DO I=1,M
       SC=1.0_SP/COSH(B*(I-Xwavemaker/DX-1)*DX)
       Eta(I,J)=A1*SC*SC+A2*SC*SC*SC*SC
       U(I,J)=A*SC*SC
       V(I,J)=ZERO
     ENDDO
     ENDDO


END SUBROUTINE INITIAL_SOLITARY_WAVE



! --------------------------------------------------
!    This is subroutine to provide solitary wave solution of Nogu equations
!    called by
!       - INITIAL_SOLITARY_WAVE
!    Last Update: 10/05/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE SUB_SLTRY(a0, h1, alp, r1, bue, ae1, ae2, au, C_ph)
    USE PARAM
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: a0,h1,alp
    REAL(SP), INTENT(OUT) :: r1,bue,ae1,ae2,au,C_ph
    REAL(SP) :: alp2,p,q,r,x,fx, fpx,rx,cph,eps,zr,g
    INTEGER :: ite

!
!--------coefficients for third order polynomial equations
!               " x**3+p*x**2+q*x+r=0  "

         alp2 = alp + 1.0_SP/3.0_SP
         eps  = a0/h1
         g=GRAV

         p = -(alp2+2.0_SP*alp*(1.0_SP+eps))/(2.0_SP*alp)
         q = eps*alp2/alp
         r = alp2/(2.0_SP*alp)

!--------Newton-Rapson's method to solve x ( >1 )for the above equation
         ite = 0
         x=1.2_SP
1        fx  = r+x*(q+x*(p+x))
         fpx = q+x*(2.0_SP*p+3.0_SP*x)
         x = x-fx/fpx

         ite = ite+1
         if (ite.gt.10) then
            !write(*,*) 'no solitary wave solution (check eps = a0/h1)'
            stop
         endif
         if (abs(fx).ge.1e-5) goto 1
         rx = sqrt(x)
         cph = sqrt(g*h1)
         r1 = rx*cph
         C_ph=rx*cph

         !write (*,*) rx, r1


!---------coefficients for solitary solutions :
!         "   u =  au/(cosh(bue*xi))**2    "
!         "  et = ae1/(cosh(bue*xi)**2+ae2/(cosh(bue*xi)**4   "

          au =  (x-1.0_SP)/(eps*rx)*cph*eps
         bue =  sqrt((x-1.0_SP)/(4.0_SP*(alp2-alp*x)))/h1
         ae1 =  (x-1.0_SP)/(eps*3.0_SP*(alp2-alp*x))*a0
         ae2 = -(x-1.0_SP)/(2.0_SP*eps)*(x-1.0_SP)*(2.0_SP*alp*x+alp2) &
                                      /(x*(alp2-alp*x))*a0
         zr  = bue*rx*cph

END SUBROUTINE SUB_SLTRY

