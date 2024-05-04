SUBROUTINE GET_COUPLING_NEEDS (COMPDA_C,KGRPNT_C,XCGRID_C,YCGRID_C,  &
                 MXK,MYK,VOQR_C,VOQ_C)
      USE OUTP_DATA                                                       
      USE SWAN_COMMON
      IMPLICIT NONE
      REAL COMPDA_C(MCGRD,MCMVAR),XCGRID_C(MXC,MYC),YCGRID_C(MXC,MYC)
      integer KGRPNT_C(MXC,MYC),MXK,MYK
      integer VOQR_C(*)
      real VOQ_C(MXK*MYK,*)

       write(*,*) MXC,MYC,MXK,MYK

END SUBROUTINE GET_COUPLING_NEEDS

SUBROUTINE SHORECIRC2SWAN
!      USE SWAN_COMMON, ONLY : IX,IY,INODE,MXC,MYC,FIRST_CALL_PASS_SC
      USE SWAN_COMMON
      USE PARAM
      USE GLOBAL, ONLY : Mloc,Nloc,U,V,ETA,Mglob,Nglob,TIME,MASK
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mglob,Nglob)::VaGlob
      REAL(SP),DIMENSION(Mloc,Nloc)::ETA_MASKED
      REAL(SP),DIMENSION(MXC,MYC)::EtaSW,USW,VSW

        IF(FIRST_CALL_PASS_SC)THEN
          FIRST_CALL_PASS_SC = .FALSE.
          IF(.NOT.ALLOCATED(GRID_ANGLE_SW2SC)) &
             ALLOCATE(GRID_ANGLE_SW2SC(MXC,MYC))
          GRID_ANGLE_SW2SC = ZERO
          DO IY=1,MYC
          DO IX=1,MXC-1
           GRID_ANGLE_SW2SC(IX,IY)= &
           ATAN2(YCGRID(IX+1,IY)-YCGRID(IX,IY),XCGRID(IX+1,IY)-XCGRID(IX,IY))
          ENDDO
          ENDDO
          DO IY=1,MYC
           GRID_ANGLE_SW2SC(MXC,IY)=GRID_ANGLE_SW2SC(MXC-1,IY)
          ENDDO

          VaGlob = ZERO
          EtaSW = ZERO
          USW = ZERO
          VSW = ZERO
        ENDIF

       IF(TIME.GE.WC_LAG)THEN 

! for serial code, some problem here because sizes of 
! ETA, MASK, EtaSW and MASK_WC_INTERACT are inconsistant 02/01/2012


        EtaSW = ETA*MASK  ! I removed *MASK_WC_INTERACT
        USW = U*MASK_WC_INTERACT
        VSW = V*MASK_WC_INTERACT



! end parallel

! assign to previous time level
        DO INDX = 2, MCGRD
          COMPDA(INDX,JWLV1)=COMPDA(INDX,JWLV2)
          COMPDA(INDX,JVX1)=COMPDA(INDX,JVX2)
          COMPDA(INDX,JVY1)=COMPDA(INDX,JVY2)
        ENDDO
         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          IF(INDX.GT.1)THEN
           COMPDA(INDX,JWLV2)=EtaSW(IX,IY)
           COMPDA(INDX,JVX2)=USW(IX,IY)*COS(GRID_ANGLE_SW2SC(IX,IY)) &
                            +VSW(IX,IY)*SIN(GRID_ANGLE_SW2SC(IX,IY))
           COMPDA(INDX,JVY2)=VSW(IX,IY)*COS(GRID_ANGLE_SW2SC(IX,IY)) &
                            -USW(IX,IY)*SIN(GRID_ANGLE_SW2SC(IX,IY))
          ENDIF
         ENDDO
         ENDDO

       ENDIF ! WC_LAG

!      IF(INODE.EQ.1)THEN
!        OPEN(2,FILE='tmp1.txt')
!         DO IY=1,MYC
!          WRITE(2,1221)(EtaSW(IX,IY),IX=1,MXC)
!         ENDDO
!        CLOSE(2)
!      ENDIF
1221     format(500f12.6)

END SUBROUTINE SHORECIRC2SWAN



SUBROUTINE SWAN2SHORECIRC
      USE SWAN_COMMON
      USE OUTP_DATA
      USE GLOBAL, ONLY : Mloc,Nloc,TIME,H,MinDepthFrc,tmp4preview
      IMPLICIT NONE

        IF(FIRST_CALL_PASS_SW)THEN
         FIRST_CALL_PASS_SW = .FALSE.
! SHORECIRC INITIALIZATION

        IF(.NOT.ALLOCATED(WaveHeightSW)) ALLOCATE(WaveHeightSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveHeightSC)) ALLOCATE(WaveHeightSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(PeakPeriodSW)) ALLOCATE(PeakPeriodSW(MXC,MYC))
        IF(.NOT.ALLOCATED(PeakPeriodSC)) ALLOCATE(PeakPeriodSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(AvePeriodSW)) ALLOCATE(AvePeriodSW(MXC,MYC))
        IF(.NOT.ALLOCATED(AvePeriodSC)) ALLOCATE(AvePeriodSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveAngleSW)) ALLOCATE(WaveAngleSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveAngleSC)) ALLOCATE(WaveAngleSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(PeakAngleSW)) ALLOCATE(PeakAngleSW(MXC,MYC))
        IF(.NOT.ALLOCATED(PeakAngleSC)) ALLOCATE(PeakAngleSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveUbottSW)) ALLOCATE(WaveUbottSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveUbottSC)) ALLOCATE(WaveUbottSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFxSW)) ALLOCATE(WaveFxSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFxSC)) ALLOCATE(WaveFxSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFySW)) ALLOCATE(WaveFySW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFySC)) ALLOCATE(WaveFySC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFluxXSW)) ALLOCATE(WaveFluxXSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFluxXSC)) ALLOCATE(WaveFluxXSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveFluxYSW)) ALLOCATE(WaveFluxYSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveFluxYSC)) ALLOCATE(WaveFluxYSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveDissSW)) ALLOCATE(WaveDissSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveDissSC)) ALLOCATE(WaveDissSC(Mloc,Nloc))

        IF(.NOT.ALLOCATED(WaveBrFraSW)) ALLOCATE(WaveBrFraSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveBrFraSC)) ALLOCATE(WaveBrFraSC(Mloc,Nloc))

        WaveHeightSW = ZERO
        WaveHeightSC = ZERO
        PeakPeriodSW = ZERO
        PeakPeriodSC = ZERO
        WaveAngleSW  = ZERO
        WaveAngleSC  = ZERO
        PeakAngleSW  = ZERO
        PeakAngleSC  = ZERO
        WaveUbottSW  = ZERO
        WaveUbottSC  = ZERO
        WaveFxSW = ZERO
        WaveFxSC = ZERO
        WaveFySW = ZERO
        WaveFySC = ZERO
        WaveFluxXSW = ZERO
        WaveFluxYSW = ZERO
        WaveFluxXSC = ZERO
        WaveFluxYSC = ZERO
        WaveDissSW = ZERO
        WaveDissSC = ZERO
        WaveBrFraSW = ZERO
        WaveBrFraSC = ZERO



        ENDIF ! end first_call_pass

! wave height and bottom orbital velocity, dissipation
! note: for dissipation in SWAN, 
! JDISS - total dissipation
! JDSXB - bottom friction dissipation
! JDSXS - surfzone breaking
! JDSXW - whitecapping dissipation
! JQB   - wave breaking fractoin
         DO IX=1,MXC
         DO IY=1,MYC
          INDX=KGRPNT(IX,IY)
          WaveHeightSW(IX,IY)=COMPDA(INDX,JHS)
          WaveUbottSW(IX,IY) =COMPDA(INDX,JUBOT)
          WaveDissSW(IX,IY) = COMPDA(INDX,JDSXS)
          WaveBrFraSW(IX,IY) = COMPDA(INDX,JQB)
         ENDDO
         ENDDO


! serial code

! wave height

        CALL DISTRIBUTE2SHORECIRC(WaveHeightSW,WaveHeightSC)

! bottom u

        CALL DISTRIBUTE2SHORECIRC(WaveUbottSW,WaveUbottSC)

! wave angle

        CALL DISTRIBUTE2SHORECIRC(WaveAngleSW,WaveAngleSC)

! wave fluxX

        CALL DISTRIBUTE2SHORECIRC(WaveFluxXSW,WaveFluxXSC)

! wave fluxY

        CALL DISTRIBUTE2SHORECIRC(WaveFluxYSW,WaveFluxYSC)

! mean period


        CALL DISTRIBUTE2SHORECIRC(AvePeriodSW,AvePeriodSC)

! peakperiod

        CALL DISTRIBUTE2SHORECIRC(PeakPeriodSW,PeakPeriodSC)

! wave dissi

        CALL DISTRIBUTE2SHORECIRC(WaveDissSW,WaveDissSC)

! wave breaking fraction

        CALL DISTRIBUTE2SHORECIRC(WaveBrFraSW,WaveBrFraSC)


! end parallel


! wave force
     IF(SHORECIRC_RUN)THEN
       IF(TIME.GE.WC_LAG)THEN 


! serial code


        WaveFxSW = WaveFxSW*MASK_WC_INTERACT


        WaveFySW = WaveFySW*MASK_WC_INTERACT

        CALL DISTRIBUTE2SHORECIRC(WaveFxSW,WaveFxSC)
        CALL DISTRIBUTE2SHORECIRC(WaveFySW,WaveFySC)
tmp4preview=WaveFySW

       ENDIF
     ENDIF ! end shorecirc_run = .true.

END SUBROUTINE SWAN2SHORECIRC



SUBROUTINE DISTRIBUTE2SHORECIRC(PHIGLOB_IN,PHI)
     USE PARAM
     USE GLOBAL, ONLY : Mloc,Nloc,Nghost,Ibeg,Iend,Jbeg,Jend,&
                      PERIODIC_X,PERIODIC_Y,Num_Transit,Mglob,Nglob
     IMPLICIT NONE

     REAL(SP),DIMENSION(Mglob,Nglob),INTENT(IN) :: PHIGLOB_IN
     REAL(SP),DIMENSION(Mglob,Nglob) :: PHIGLOB
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: PHI
     REAL(SP) :: TMP_L, TMP_R

     PHIGLOB=PHIGLOB_IN

     IF(PERIODIC_X)THEN
       DO J=1,NGlob
        TMP_L=PHIGLOB(MGlob-Num_Transit+1,J)
        TMP_R=PHIGLOB(Num_Transit,J)
        DO I=1,Num_Transit
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(I+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO I=MGlob-Num_Transit,MGlob
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(I-MGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

     IF(PERIODIC_Y)THEN
       DO I=1,MGlob
        TMP_L=PHIGLOB(I,NGlob-Num_Transit+1)
        TMP_R=PHIGLOB(I, Num_Transit)
        DO J=1,Num_Transit
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(J+Num_Transit-1)/(2.0*Num_Transit)
        ENDDO
        DO J=NGlob-Num_Transit,NGlob
          PHIGLOB(I,J)=TMP_L+(TMP_R-TMP_L)*(J-NGlob+Num_Transit)/(2.0*Num_Transit)
        ENDDO
       ENDDO
     ENDIF

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       PHI(I,J)=PHIGLOB(I-Ibeg+1,J-Jbeg+1)
     ENDDO
     ENDDO

! ghost cells
     DO I=Ibeg,Iend
       DO J=1,Nghost
        PHI(I,J)=PHI(I,Jbeg)
       ENDDO
       DO J=Jend+1,Nloc
        PHI(I,J)=PHI(I,Jend)
       ENDDO
     ENDDO

     DO J=1,Nloc
       DO I=1,Nghost
        PHI(I,J)=PHI(Ibeg,J)
       ENDDO
       DO I=Iend+1,Mloc
        PHI(I,J)=PHI(Iend,J)
       ENDDO
     ENDDO 


END SUBROUTINE DISTRIBUTE2SHORECIRC








