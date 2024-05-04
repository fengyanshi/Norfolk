! ----------------------------------------------------
!    This is subroutine to read input.txt
!  called by 
!        MAIN
!    Last Update: 05/07/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------

SUBROUTINE READ_INPUT
    USE GLOBAL
    USE Input_Util
    USE PASS






    IMPLICIT NONE
    CHARACTER(LEN=80) FILE_NAME
    INTEGER::LINE,I_tmp
    INTEGER :: ierr
    CHARACTER(LEN=80) :: WHAT





      OPEN(3,FILE='LOG.txt')   

! read everything from input.txt







      FILE_NAME='INPUT'


! title
      CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
      IF(ierr==1)THEN
        !write(*,*) 'No TITLE in ', FILE_NAME, 'use default'
        TITLE='---TEST RUN---'
      ENDIF





      WRITE(3,*)'---- LOG FILE ---'
      WRITE(3,*)TITLE
      WRITE(3,*)' --------------input start --------------'     







! wave-current interaction
      CALL GET_LOGICAL_VAL(SWAN_RUN,FILE_NAME,'SWAN_RUN',line)
      CALL GET_LOGICAL_VAL(SHORECIRC_RUN,FILE_NAME,'SHORECIRC_RUN',line)

      CALL GET_INTEGER_VAL(WC_BOUND_WEST,FILE_NAME,'WC_BOUND_WEST',line)
      CALL GET_INTEGER_VAL(WC_BOUND_EAST,FILE_NAME,'WC_BOUND_EAST',line)
      CALL GET_INTEGER_VAL(WC_BOUND_SOUTH,FILE_NAME,'WC_BOUND_SOUTH',line)
      CALL GET_INTEGER_VAL(WC_BOUND_NORTH,FILE_NAME,'WC_BOUND_NORTH',line)
      CALL GET_Float_VAL(WC_LAG,FILE_NAME,'WC_LAG',line)  







    







      WRITE(3,'(A16,I4)') 'WC_BOUND_WEST=',WC_BOUND_WEST
      WRITE(3,'(A16,I4)') 'WC_BOUND_EAST=',WC_BOUND_EAST
      WRITE(3,'(A16,I4)') 'WC_BOUND_SOUTH=',WC_BOUND_SOUTH
      WRITE(3,'(A16,I4)') 'WC_BOUND_NORTH=',WC_BOUND_NORTH


! dimension
      CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
      CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)



      WRITE(3,'(A7,I3,A7,I3)') 'Mglob=',Mglob,'Nglob=', Nglob


! grid 


        CALL GET_STRING_VAL(X_FILE,FILE_NAME,'X_FILE',line,ierr)
        CALL GET_STRING_VAL(Y_FILE,FILE_NAME,'Y_FILE',line,ierr)




        WRITE(3,'(A9,A50)') 'X_FILE= ', X_FILE
        WRITE(3,'(A9,A50)') 'Y_FILE= ', Y_FILE

    


! coriolis
      CALL GET_LOGICAL_VAL(CORI_CONSTANT,FILE_NAME,'CORI_CONSTANT',line)   
     if (CORI_CONSTANT) then
      CALL GET_Float_VAL(LATITUDE,FILE_NAME,'LATITUDE',line)
     else
        CALL GET_STRING_VAL(LATITUDE_FILE,FILE_NAME,'LATITUDE_FILE',line,ierr) 
     endif

! result folder
      CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)



      WRITE(3,'(A15,A50)')'RESULT_FOLDER:', RESULT_FOLDER

! station files
      CALL GET_INTEGER_VAL(NumberStations,FILE_NAME,'NumberStations',line)
      IF(NumberStations>0)THEN
      CALL GET_STRING_VAL(STATIONS_FILE,FILE_NAME,'STATIONS_FILE',line,ierr)
      ENDIF
! depth 
      CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)



      WRITE(3,'(A12,A50)')'DEPTH_TYPE:', DEPTH_TYPE

      IF(DEPTH_TYPE(1:3)=='DAT')THEN
        CALL GET_STRING_VAL(DEPTH_FILE,FILE_NAME,'DEPTH_FILE',line,ierr)



      WRITE(3,'(A12,A50)')'DEPTH_FILE:', DEPTH_FILE

      ENDIF
      IF(DEPTH_TYPE(1:3)=='FLA')THEN
      CALL GET_Float_VAL(DEPTH_FLAT,FILE_NAME,'DEPTH_FLAT',line) 



      WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT 

      ENDIF
      IF(DEPTH_TYPE(1:3)=='SLO')THEN
      CALL GET_Float_VAL(DEPTH_FLAT,FILE_NAME,'DEPTH_FLAT',line) 
      CALL GET_Float_VAL(SLP,FILE_NAME,'SLP',line) 
      CALL GET_Float_VAL(Xslp,FILE_NAME,'Xslp',line) 





      WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT 
      WRITE(3,'(A5,F12.2)')'SLP=', SLP
      WRITE(3,'(A6,F12.2)')'Xslp=', Xslp  

      ENDIF
! time
      CALL GET_Float_VAL(PLOT_INTV,FILE_NAME,'PLOT_INTV',line)
      CALL GET_Float_VAL(PLOT_INTV_STATION,FILE_NAME,'PLOT_INTV_STATION',line)
      CALL GET_Float_VAL(SCREEN_INTV,FILE_NAME,'SCREEN_INTV',line)




      WRITE(3,'(A12,F12.2)')'TOTAL_TIME=', TOTAL_TIME
      WRITE(3,'(A12,F12.2)')'PLOT_INTV= ', PLOT_INTV
      WRITE(3,'(A13,F12.2)')'SCREEN_INTV=', SCREEN_INTV
      WRITE(3,'(A15,F12.2)')'HOTSTART_INTV=', HOTSTART_INTV


! stationary mode

! end stationary

! initial uvz
      CALL GET_LOGICAL_VAL(INI_UVZ,FILE_NAME,'INI_UVZ',line)
      IF(INI_UVZ)THEN
        CALL GET_STRING_VAL(ETA_FILE,FILE_NAME,'ETA_FILE',line,ierr)
        CALL GET_STRING_VAL(U_FILE,FILE_NAME,'U_FILE',line,ierr)
        CALL GET_STRING_VAL(V_FILE,FILE_NAME,'V_FILE',line,ierr)
       ENDIF
! open boundary conditions
      CALL GET_LOGICAL_VAL(ETA_CLAMPED,FILE_NAME,'ETA_CLAMPED',line)
      IF(ETA_CLAMPED)THEN
        CALL GET_STRING_VAL(TIDE_FILE,FILE_NAME,'TIDE_FILE',line,ierr)
      ENDIF



        WRITE(3,'(A12,A50)') 'TIDE_FILE= ', TIDE_FILE

! tidal flux
      CALL GET_LOGICAL_VAL(FLUX_TIDE,FILE_NAME,'FLUX_TIDE',line)
      IF(FLUX_TIDE)THEN
        CALL GET_STRING_VAL(FLUX_TIDE_FILE,FILE_NAME,'FLUX_TIDE_FILE',line,ierr)
      ENDIF



        WRITE(3,'(A16,A50)') 'FLUX_TIDE_FILE=', FLUX_TIDE_FILE



      CALL GET_LOGICAL_VAL(FLUX_CLAMPED,FILE_NAME,'FLUX_CLAMPED',line)
      IF(FLUX_CLAMPED)THEN
        CALL GET_STRING_VAL(FLUX_FILE,FILE_NAME,'FLUX_FILE',line,ierr)
      ENDIF



        WRITE(3,'(A12,A50)') 'FLUX_FILE= ', FLUX_FILE

! flux
      CALL GET_LOGICAL_VAL(FLUX_CLAMPED,FILE_NAME,'FLUX_CLAMPED',line)
      IF(FLUX_CLAMPED)THEN
        CALL GET_STRING_VAL(FLUX_FILE,FILE_NAME,'FLUX_FILE',line,ierr)
      ENDIF



        WRITE(3,'(A12,A50)') 'FLUX_FILE= ', FLUX_FILE


      CALL GET_LOGICAL_VAL(WindForce,FILE_NAME,'WindForce',line)
      IF(WindForce)THEN
        CALL GET_STRING_VAL(WIND_FILE,FILE_NAME,'WIND_FILE',line,ierr)
        CALL GET_Float_VAL(Cdw,FILE_NAME,'Cdw',line)
      ENDIF



      WRITE(3,'(A7,A50)')'Cdw = ', Cdw


! wavemaker
      CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)



      WRITE(3,'(A11,A50)')'WAVEMAKER:', WAVEMAKER

        IF(WaveMaker(1:7)=='LEF_SOL')THEN
          CALL GET_Float_VAL(AMP_SOLI,FILE_NAME,'AMP',line)
          CALL GET_Float_VAL(DEP_SOLI,FILE_NAME,'DEP',line)
          CALL GET_Float_VAL(LAG_SOLI,FILE_NAME,'LAGTIME',line)





      WRITE(3,'(A10,F12.2)')'AMP_SOLI=', AMP_SOLI
      WRITE(3,'(A10,F12.2)')'DEP_SOLI=', DEP_SOLI
      WRITE(3,'(A10,F12.2)')'LAG_SOLI=', LAG_SOLI

        ENDIF

        IF(WaveMaker(1:7)=='INI_SOL')THEN
          CALL GET_Float_VAL(AMP_SOLI,FILE_NAME,'AMP',line)
          CALL GET_Float_VAL(DEP_SOLI,FILE_NAME,'DEP',line)
          CALL GET_Float_VAL(XWAVEMAKER,FILE_NAME,'XWAVEMAKER',line)




      WRITE(3,'(A10,F12.2)')'AMP_SOLI=', AMP_SOLI
      WRITE(3,'(A10,F12.2)')'DEP_SOLI=', DEP_SOLI

        ENDIF
        IF(WaveMaker(1:6)=='N_WAVE')THEN
          CALL GET_Float_VAL(x1_Nwave,FILE_NAME,'x1_Nwave',line)
          CALL GET_Float_VAL(x2_Nwave,FILE_NAME,'x2_Nwave',line)
          CALL GET_Float_VAL(a0_Nwave,FILE_NAME,'a0_Nwave',line)
          CALL GET_Float_VAL(gamma_Nwave,FILE_NAME,'gamma_Nwave',line)
          CALL GET_Float_VAL(dep_Nwave,FILE_NAME,'dep_Nwave',line)







      WRITE(3,'(A10,F12.2)')'x1_Nwave=', x1_Nwave
      WRITE(3,'(A10,F12.2)')'x2_Nwave=', x2_Nwave
      WRITE(3,'(A10,F12.2)')'a0_Nwave=', a0_Nwave
      WRITE(3,'(A13,F12.2)')'gamma_Nwave=', gamma_Nwave
      WRITE(3,'(A11,F12.2)')'dep_Nwave=', dep_Nwave

        ENDIF

        IF(WaveMaker(1:7)=='INI_REC')THEN
          CALL GET_Float_VAL(AMP_SOLI,FILE_NAME,'AMP',line)
          CALL GET_Float_VAL(Xc,FILE_NAME,'Xc',line)
          CALL GET_Float_VAL(Yc,FILE_NAME,'Yc',line)
          CALL GET_Float_VAL(WID,FILE_NAME,'WID',line)






      WRITE(3,'(A10,F12.2)')'AMP     =', AMP_SOLI
      WRITE(3,'(A10,F12.2)')'Xc      =', Xc
      WRITE(3,'(A10,F12.2)')'Yc      =', Yc
      WRITE(3,'(A10,F12.2)')'WID     =', WID

        ENDIF

! periodic boundary condition
      CALL GET_LOGICAL_VAL(PERIODIC_X,FILE_NAME,'PERIODIC_X',line)
      CALL GET_LOGICAL_VAL(PERIODIC_Y,FILE_NAME,'PERIODIC_Y',line)
      CALL GET_INTEGER_VAL(Num_Transit,FILE_NAME,'Num_Transit',line)






      WRITE(3,'(A13,A50)')'PERIODIC_X:', PERIODIC_X
      WRITE(3,'(A13,A50)')'PERIODIC_Y:', PERIODIC_Y
      WRITE(3,'(A15,I4)')'Num_Transit:', Num_Transit


! sponge layer
      CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line)



      WRITE(3,'(A11,A50)')'SPONGE_ON:', SPONGE_ON

      IF(SPONGE_ON)THEN
        CALL GET_Float_VAL(Sponge_west_width,FILE_NAME,'Sponge_west_width',line)
        CALL GET_Float_VAL(Sponge_east_width,FILE_NAME,'Sponge_east_width',line)
        CALL GET_Float_VAL(Sponge_south_width,FILE_NAME,'Sponge_south_width',line)
        CALL GET_Float_VAL(Sponge_north_width,FILE_NAME,'Sponge_north_width',line)
        CALL GET_Float_VAL(R_sponge,FILE_NAME,'R_sponge',line)
        CALL GET_Float_VAL(A_sponge,FILE_NAME,'A_sponge',line)








        WRITE(3,'(A20,F12.2)')'Sponge_west_width =', Sponge_west_width
        WRITE(3,'(A20,F12.2)')'Sponge_east_width =', Sponge_east_width
        WRITE(3,'(A20,F12.2)')'Sponge_south_width=', Sponge_south_width
        WRITE(3,'(A20,F12.2)')'Sponge_north_width=', Sponge_north_width
        WRITE(3,'(A20,F12.2)')'R_sponge          =', R_sponge
        WRITE(3,'(A20,F12.2)')'A_sponge          =', A_sponge

      ENDIF

! obstacle structures
      CALL GET_STRING_VAL(OBSTACLE_FILE,FILE_NAME,'OBSTACLE_FILE',line,ierr)
      IF(ierr==1)THEN
        OBSTACLE=.FALSE.



      WRITE(3,'(A15,A5)')'OBSTACLE_FILE:', 'NO'

      ELSE
        OBSTACLE=.TRUE.



      WRITE(3,'(A15,A50)')'OBSTACLE_FILE:', OBSTACLE_FILE

      ENDIF

! physics



      WRITE(3,'(A8)')'Physics'

      CALL GET_Float_VAL(Cd,FILE_NAME,'Cd',line)
      CALL GET_Float_VAL(Manning,FILE_NAME,'Manning',line)
      CALL GET_Float_VAL(nu_bkgd,FILE_NAME,'nu_bkgd',line)

      CALL GET_LOGICAL_VAL(FRC_MANNING_DATA,FILE_NAME,'FRC_MANNING_DATA',line)
      IF(FRC_MANNING_DATA)THEN
        CALL GET_STRING_VAL(FRC_FILE,FILE_NAME,'FRC_FILE',line,ierr)
      ENDIF






       WRITE(3,'(A13,F12.2)')'Cd         =', Cd
       WRITE(3,'(A13,F12.2)')'Manning    =', Manning
       WRITE(3,'(A13,F12.2)')'nu_bkgd    =', nu_bkgd

! numerics schemes
      CALL GET_STRING_VAL(Time_Scheme,FILE_NAME,'Time_Scheme',line,ierr)
      IF(ierr==1)THEN
        !write(*,*) 'Please define Time_Scheme in ', FILE_NAME



      WRITE(3,'(A13,A50)')'TIME_SCHEME:', 'NOT DEFINED, STOP'

        STOP
      ENDIF



      WRITE(3,'(A13,A50)')'TIME_SCHEME:', TIME_SCHEME

      CALL GET_STRING_VAL(CONSTR,FILE_NAME,'CONSTRUCTION',line,ierr)
      IF(ierr==1)THEN
        !write(*,*) 'No definition of CONSTRUCTION in ', FILE_NAME, 'use default'



      WRITE(3,'(A14,A50)')'CONSTRUCTION', 'NOT DEFINED, USE DEFAULT'

        CONSTR='HLLC'
      ENDIF



      WRITE(3,'(A14,A50)')'CONSTRUCTION:', CONSTR

      CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
      IF(ierr==1)THEN
        !write(*,*) 'No definition of HIGH_ORDER in ', FILE_NAME, 'use default'



      WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'

        HIGH_ORDER='FOURTH'        
      ENDIF



      WRITE(3,'(A12,A50)')'HIGH_ORDER:', HIGH_ORDER

! CFL
      CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)



      WRITE(3,'(A5,F12.2)')'CFL=', CFL

! Froude Number Cap
      CALL GET_Float_VAL(FroudeCap,FILE_NAME,'FroudeCap',line)



      WRITE(3,'(A5,F12.2)')'FroudeCap=', FroudeCap

! MinDepth etc
      CALL GET_Float_VAL(MinDepth,FILE_NAME,'MinDepth',line)



      WRITE(3,'(A10,F12.6)')'MinDepth=', MinDepth

      CALL GET_Float_VAL(MinDepthFrc,FILE_NAME,'MinDepthFrc',line)



      WRITE(3,'(A13,F12.2)')'MinDepthFrc=', MinDepthFrc














! end of sediment


      CALL GET_STRING_VAL(COUPLING_FILE,FILE_NAME,'COUPLING_FILE',line,ierr)



      WRITE(3,'(A15,A50)')'COUPLING_FILE:', COUPLING_FILE



! output parameters
      CALL GET_LOGICAL_VAL(OUT_DEPTH,FILE_NAME,'DEPTH_OUT',line)
      CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'U',line)
      CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'V',line)
      CALL GET_LOGICAL_VAL(OUT_ETA,FILE_NAME,'ETA',line)
      CALL GET_LOGICAL_VAL(OUT_HS,FILE_NAME,'HS',line)
      CALL GET_LOGICAL_VAL(OUT_MASK,FILE_NAME,'MASK',line)
      CALL GET_LOGICAL_VAL(OUT_SourceX,FILE_NAME,'SourceX',line)
      CALL GET_LOGICAL_VAL(OUT_SourceY,FILE_NAME,'SourceY',line)
      CALL GET_LOGICAL_VAL(OUT_TMP,FILE_NAME,'TMP',line)
      CALL GET_LOGICAL_VAL(OUT_Per,FILE_NAME,'PER',line)
      CALL GET_LOGICAL_VAL(OUT_WDIR,FILE_NAME,'WDIR',line)
      CALL GET_LOGICAL_VAL(OUT_Wdis,FILE_NAME,'Wdis',line)
      CALL GET_LOGICAL_VAL(OUT_WBV,FILE_NAME,'WBV',line)
      CALL GET_LOGICAL_VAL(OUT_WFC,FILE_NAME,'WFC',line)
      CALL GET_LOGICAL_VAL(OUT_3D,FILE_NAME,'UV3D',line)
      CALL GET_LOGICAL_VAL(OUT_Qw,FILE_NAME,'Qstk',line)


! 



      WRITE(3,*)' --------------input end --------------' 


END SUBROUTINE READ_INPUT

! -------------
!    Writes station data
! Fengyan Shi modified based on Jeff Harris' for Spherical
! here simply specify grid number i and j instead of x and y
! 04/14/2011
! -------------
SUBROUTINE STATIONS
     USE GLOBAL
     USE PASS
     IMPLICIT NONE

     INTEGER :: iunit
     REAL(SP) :: dum1,dum2
     CHARACTER(LEN=80)::FILE_NAME=''
     CHARACTER(LEN=80)::TMP_NAME=''
     CHARACTER(LEN=80)::FDIR=''

! initialize stations
     FDIR=TRIM(RESULT_FOLDER)
     if (icount.eq.0) then
       ALLOCATE(ista(NumberStations),&
                jsta(NumberStations),&
                nsta(NumberStations))
! calculate how many output components
              
       open(100,FILE=TRIM(STATIONS_FILE))
       do i=1,NumberStations
          read(100,*) dum1,dum2

          ista(i) = Nghost+dum1
          jsta(i) = Nghost+dum2
          if ((ista(i).ge.Ibeg).and.(ista(i).le.Iend).and.&
              (jsta(i).ge.Jbeg).and.(jsta(i).le.Jend)) then
             nsta(i) = 1
             write(file_name(1:4),'(I4.4)') i
             TMP_NAME = TRIM(FDIR)//'sta_'//TRIM(FILE_NAME)
             iunit=100+i
             open(iunit,FILE=TMP_NAME)
          else
             nsta(i) = 0
          endif

       enddo
     endif

! write to stations

     do i=1,NumberStations
       if (nsta(i).eq.1) then
          iunit=100+i
          write (iunit,'(20E20.6)') time, eta(ista(i),jsta(i)),&
                          u(ista(i),jsta(i)),v(ista(i),jsta(i)),&
                          depth(ista(i),jsta(i))
       endif
     enddo

! close station files
     if (TIME.ge.TOTAL_TIME) then
       do i=1,NumberStations
          if (nsta(i).eq.1) then
             iunit=100+i
             close(iunit)
          endif
       enddo
     endif

END SUBROUTINE STATIONS

! ----------------------------------------------------
!    This is subroutine for preview
!  called by 
!        MAIN
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE PREVIEW
     USE GLOBAL
     USE PASS



     IMPLICIT NONE

     CHARACTER(LEN=80)::FILE_NAME=''
     CHARACTER(LEN=80)::TMP_NAME=''
     CHARACTER(LEN=80)::FDIR=''
     INTEGER :: numprint,VTYPE

     FDIR=TRIM(RESULT_FOLDER)

     ICOUNT=ICOUNT+1







        WRITE(*,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time

102     FORMAT(A20,I4,A14,F12.3,A2,F12.3)





        numprint=icount          !    +myid*1000
        itmp1=mod(numprint/1000,10)
        itmp2=mod(numprint/100,10)
        itmp3=mod(numprint/10,10)
        itmp4=mod(numprint,10)


        write(file_name(1:1),'(I1)')itmp1
        write(file_name(2:2),'(I1)')itmp2
        write(file_name(3:3),'(I1)')itmp3
        write(file_name(4:4),'(I1)')itmp4


     IF(OUT_DEPTH)THEN
        TMP_NAME = TRIM(FDIR)//'dep_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,DEPTH)
     ENDIF


   IF(SHORECIRC_RUN)THEN

     IF(OUT_ETA)THEN
        TMP_NAME = TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Eta)
     ENDIF

     IF(OUT_U)THEN
        TMP_NAME = TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,U)
        TMP_NAME = TRIM(FDIR)//'us_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Usurf)
        TMP_NAME = TRIM(FDIR)//'ub_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Ubott)

        TMP_NAME = TRIM(FDIR)//'un_'//TRIM(FILE_NAME)
        Int2Flo=0.5_SP*(UnL(1:Mloc,1:Nloc)+UnR(1:Mloc,1:Nloc))
        call PutFile(TMP_NAME,Int2flo)

     ENDIF

     IF(OUT_V)THEN
        TMP_NAME = TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,V)
        TMP_NAME = TRIM(FDIR)//'vs_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Vsurf)
        TMP_NAME = TRIM(FDIR)//'vb_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Vbott)

        Int2Flo=0.5_SP*(VnL(1:Mloc,1:Nloc)+VnR(1:Mloc,1:Nloc))
        TMP_NAME = TRIM(FDIR)//'vn_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,Int2Flo)

     ENDIF

     IF(OUT_MASK)THEN
        TMP_NAME = TRIM(FDIR)//'mask_'//TRIM(FILE_NAME)
        Int2Flo=MASK
        call PutFile(TMP_NAME,Int2Flo)
     ENDIF



     IF(OUT_3D)THEN
        TMP_NAME = TRIM(FDIR)//'f1x_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,ff11)
        TMP_NAME = TRIM(FDIR)//'f1y_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,ff12)
        TMP_NAME = TRIM(FDIR)//'d1x_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,dd11)
        TMP_NAME = TRIM(FDIR)//'d1y_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,dd12)
        TMP_NAME = TRIM(FDIR)//'e1x_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,ee11)
        TMP_NAME = TRIM(FDIR)//'e1y_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,ee12)       
     ENDIF

     IF(OUT_Qw)THEN
        TMP_NAME = TRIM(FDIR)//'qwx_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveFluxXSC)
        TMP_NAME = TRIM(FDIR)//'qwy_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveFluxYSC)  
     ENDIF

     IF(OUT_SourceY)THEN
        TMP_NAME = TRIM(FDIR)//'SourceY_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,SourceY)
     ENDIF
     IF(OUT_SourceX)THEN
        TMP_NAME = TRIM(FDIR)//'SourceX_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,SourceX)
     ENDIF

   ENDIF ! end shorecirc output

! output for swan !!!!!!!!
    IF(SWAN_RUN)THEN
     IF(OUT_HS)THEN
        TMP_NAME = TRIM(FDIR)//'hs_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveHeightSC)
     ENDIF

     IF(OUT_PER)THEN
        TMP_NAME = TRIM(FDIR)//'per_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,PeakPeriodSC)

        TMP_NAME = TRIM(FDIR)//'mper_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,AvePeriodSC)
     ENDIF

     IF(OUT_Wdis)THEN
        TMP_NAME = TRIM(FDIR)//'Wdis_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveDissSC)
        TMP_NAME = TRIM(FDIR)//'Wbrk_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveBrFraSC)
     ENDIF

     IF(OUT_WDIR)THEN
        TMP_NAME = TRIM(FDIR)//'wdir_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveAngleSC)
     ENDIF

     IF(OUT_WBV)THEN
        TMP_NAME = TRIM(FDIR)//'wbv_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveUbottSC)
     ENDIF

     IF(OUT_WFC)THEN
        TMP_NAME = TRIM(FDIR)//'wfx_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveFxSC)
        TMP_NAME = TRIM(FDIR)//'wfy_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveFySC)
     ENDIF
    ENDIF ! end swan output














210   FORMAT(5000I3)


     IF(OUT_TMP)THEN
        TMP_NAME = TRIM(FDIR)//'tmp_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,tmp4preview)
     ENDIF





101   continue

END SUBROUTINE PREVIEW


SUBROUTINE GetFile (FILE,PHI)
     USE PARAM
     USE GLOBAL, ONLY : Mloc,Nloc,Nghost
     IMPLICIT NONE
     CHARACTER(LEN=80) FILE
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: PHI

      OPEN(1,FILE=TRIM(FILE))
       DO J=Nghost+1,Nloc-Nghost
        READ(1,*)(PHI(I,J),I=Nghost+1,Mloc-Nghost)
       ENDDO
      CLOSE(1)

! ghost cells for initial uvz only
        DO I=Nghost+1,Mloc-Nghost
           DO J=1,Nghost
              PHI(I,J)=PHI(I,Nghost+1)
           ENDDO
           DO J=Nloc-Nghost+1,Nloc
              PHI(I,J)=PHI(I,Nloc-Nghost)
           ENDDO
        ENDDO
        DO J=1,Nloc
           DO I=1,Nghost
              PHI(I,J)=PHI(Nghost+1,J)
           ENDDO
           DO I=Mloc-Nghost+1,Mloc
              PHI(I,J)=PHI(Mloc-Nghost,J)
           ENDDO
        ENDDO

END SUBROUTINE Getfile



SUBROUTINE PutFile(FILE,PHI)
     USE PARAM
     USE GLOBAL
     IMPLICIT NONE
     REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI
     CHARACTER(LEN=80) FILE
        OPEN(1,FILE=TRIM(FILE))





        DO J=Nghost+1,Nloc-Nghost
           WRITE(1,100)(real(PHI(I,J)),I=Nghost+1,Mloc-Nghost)
        ENDDO


100  FORMAT(5000E16.6)
        CLOSE(1)
END SUBROUTINE PutFile









! --------------------------------------------------
!    This is subroutine to read hot start data and initialize other 
!    variables
!    called by
!       MAIN
!    
!    Last Update: 10/21/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------
SUBROUTINE READ_HOTSTART_DATA
     USE GLOBAL
     IMPLICIT NONE
     CHARACTER(LEN=80)::FDIR=''
     CHARACTER(LEN=80)::WHAT
     INTEGER :: Icontr
     INTEGER :: RMloc,RNloc,RNghost
     REAL(SP) :: RTOTAL_TIME,RPLOT_INTV,RSCREEN_INTV,RHOTSTART_INTV
     REAL(SP) :: RDX,RDY
     LOGICAL :: RDISPERSION,RSPONGE_ON
     REAL(SP) :: RGamma1,RSWE_ETA_DEP
     REAL(SP) :: RGamma2

     CHARACTER(LEN=80) RTime_Scheme
     CHARACTER(LEN=80) RCONSTR
     CHARACTER(LEN=80) RHIGH_ORDER
     REAL(SP) :: RCFL,RFroudeCap
     REAL(SP) :: RMinDepth,RMinDepthfrc
     REAL(SP) :: RSponge_west_width,RSponge_east_width, &
                 RSponge_south_width,RSponge_north_width, &
                 RR_sponge,RA_sponge
     CHARACTER(LEN=80)::FILE_NAME=''


        itmp1=mod(FileNumber_HOTSTART/1000,10)
        itmp2=mod(FileNumber_HOTSTART/100,10)
        itmp3=mod(FileNumber_HOTSTART/10,10)
        itmp4=mod(FileNumber_HOTSTART,10)


        !write(file_name(1:1),'(I1)')itmp1
        !write(file_name(2:2),'(I1)')itmp2
        !write(file_name(3:3),'(I1)')itmp3
        !write(file_name(4:4),'(I1)')itmp4

     FDIR=TRIM(RESULT_FOLDER)
     
     !WRITE(*,*)'READ IN HOT START DATA ...'
     
     OPEN(4,FILE=TRIM(FDIR)//'hotstart.'//TRIM(FILE_NAME))
    
     CLOSE(4)

     !WRITE(*,*)'READ HOT START DATA OVER!'
    
     CLOSE(4)

END SUBROUTINE READ_HOTSTART_DATA






SUBROUTINE WRITE_HOTSTART_DATA
     USE GLOBAL
     IMPLICIT NONE
     CHARACTER(LEN=80)::FDIR=''
     CHARACTER(LEN=80)::FILE_NAME=''


     FDIR=TRIM(RESULT_FOLDER)

     ICOUNT_HOTSTART=ICOUNT_HOTSTART+1

        WRITE(*,102)'HOTSTART FILE NO.', icount_hotstart
102     FORMAT(A20,I4)

        itmp1=mod(icount_hotstart/1000,10)
        itmp2=mod(icount_hotstart/100,10)
        itmp3=mod(icount_hotstart/10,10)
        itmp4=mod(icount_hotstart,10)


        !write(file_name(1:1),'(I1)')itmp1
        !write(file_name(2:2),'(I1)')itmp2
        !write(file_name(3:3),'(I1)')itmp3
        !write(file_name(4:4),'(I1)')itmp4

     !WRITE(*,*)'!WRITE OUT HOT START DATA ...'
     
     OPEN(4,FILE=TRIM(FDIR)//'hotstart.'//TRIM(FILE_NAME))
! -- dimension
     !WRITE(4,*)'# Mloc,Nloc,Mloc1,Nloc1,Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,Nghost'
     !WRITE(4,*)Mloc,Nloc,Mloc1,Nloc1,Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,Nghost
! -- time
     !WRITE(4,*)'# TIME,TOTAL_TIME,PLOT_INTV,PLOT_COUNT'
     !WRITE(4,*)'# SCREEN_INTV,SCREEN_COUNT,HOTSTART_INTV'
     !WRITE(4,*)TIME,TOTAL_TIME,PLOT_INTV,PLOT_COUNT,&
!               SCREEN_INTV,SCREEN_COUNT,HOTSTART_INTV
     !WRITE(4,*)'# ICOUNT'
     !WRITE(4,*)ICOUNT
     !WRITE(4,*)'# ICOUNT_HOTSTART'
     !WRITE(4,*)ICOUNT_HOTSTART
! -- physics
! DISPERSION
     !WRITE(4,*)'# DISPERSION'
     !WRITE(4,*)DISPERSION

     !WRITE(4,*)'# Gamma1,a1,a2,b1,b2,SWE_ETA_DEP'
     !WRITE(4,*)Gamma1,a1,a2,b1,b2,SWE_ETA_DEP
! -- numerics
     !WRITE(4,*)'# Time_Scheme'
     !WRITE(4,*)Time_Scheme
     
     !WRITE(4,*)'# HIGH_ORDER,CONSTR'
     !WRITE(4,*)HIGH_ORDER,CONSTR

     !WRITE(4,*)'# CFL,FroudeCap'
     !WRITE(4,*)CFL,FroudeCap

! -- wet-dry
     !WRITE(4,*)'# MinDepth,MinDepthfrc'
     !WRITE(4,*)MinDepth,MinDepthfrc

! -- depth

      !WRITE(4,*)'# DEPTH'
      !WRITE(4,*)((DEPTH(I,J),I=1,Mloc),J=1,Nloc)
      !WRITE(4,*)'# DEPTHx'
      !WRITE(4,*)((DEPTHx(I,J),I=1,Mloc),J=1,Nloc)
      !WRITE(4,*)'# DEPTHy'
      !WRITE(4,*)((DEPTHy(I,J),I=1,Mloc),J=1,Nloc)
! variables
      !WRITE(4,*)'# U'
      !WRITE(4,*)((U(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# V'
      !WRITE(4,*)((V(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# U0'
      !WRITE(4,*)((U0(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# V0'
      !WRITE(4,*)((V0(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# Ubar'
      !WRITE(4,*)((Ubar(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# Vbar'
      !WRITE(4,*)((Vbar(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# ETA'
      !WRITE(4,*)((ETA(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# H'
      !WRITE(4,*)((H(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# MASK'
      !WRITE(4,*)((MASK(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# MASK9'
      !WRITE(4,*)((MASK9(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# MASK_STRUC'
      !WRITE(4,*)((MASK_STRUC(I,J),I=1,Mloc),J=1,Nloc)

      !WRITE(4,*)'# SPONGE_ON'
      !WRITE(4,*) SPONGE_ON
  
      !WRITE(4,*)'# SPONGE Width'
      !WRITE(4,*)Sponge_west_width,Sponge_east_width, &
!                 Sponge_south_width,Sponge_north_width
  
      !WRITE(4,*)'# SPONGE R and A'
      !WRITE(4,*)R_sponge,A_sponge
      
    
     CLOSE(4)

END SUBROUTINE !WRITE_HOTSTART_DATA


