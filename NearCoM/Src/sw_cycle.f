         SUBROUTINE SWAN_CYCLE
         USE SWAN_COMMON
         IMPLICIT NONE
         CHARACTER PTYPE, PNAME *8, COMPUT *4, DTTIWR*18  


!           --- synchronize nodes

            CALL SWSYNC
                                                                          40.30
!
!           --- update boundary conditions and input fields
!
!TIMG            CALL SWTSTA(7)                                                40.23
            CALL SNEXTI ( BSPECS, BGRIDP, COMPDA, AC1   , AC2   ,         40.31
     &                    SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT,         40.31
     &                    XYTST , DEPTH , WLEVL , FRIC  , UXB   ,         40.31
     &                    UYB   , WXI   , WYI   )                         40.31
!TIMG            CALL SWTSTO(7)                                                40.23
!
!           --- synchronize nodes

            CALL SWSYNC                                                   40.30

            IF (COMPUT.NE.'NOCO' .AND. IT.GT.0) THEN                      40.00

              SAVITE = ITEST                                              30.21
              IF (ICOTES .GT. ITEST) ITEST = ICOTES
!
!             --- compute action density for current time step
!
!TIMG              CALL SWTSTA(8)                                              40.23
              CALL SWCOMP( AC1   , AC2   , COMPDA, SPCDIR, SPCSIG,        40.31
     &                     XYTST , IT    , KGRPNT, XCGRID, YCGRID,        40.31
     &                     CROSS )                                        40.31
!TIMG              CALL SWTSTO(8)                                              40.23
!
!             --- set ICOND=4 for stationary computation, for next
!                 (stationary) COMPUTE command                            40.13
              ICOND = 4                                                   40.13
!
!             --- check whether computed significant wave height at       32.01
!                 boundary differs from prescribed value given in         32.01
!                 boundary command values of incident Hs                  32.01
!
              IF ( BNDCHK ) THEN                                          32.01
                CALL HSOBND ( AC2, SPCSIG, COMPDA(1,JHSIBC), KGRPNT )     40.31
              ENDIF                                                       32.01
!
              ITEST = SAVITE                                              30.21

            ENDIF
            
!
            IF ( IT.EQ.IT0 .AND. .NOT.ALLOCATED(OURQT) ) THEN             40.51 40.31
               ALLOCATE (OURQT(MAX_OUTP_REQ))                             40.51 40.30
               OURQT = -9999.                                             40.51 40.30
            ENDIF                                                         40.00
!
            SAVITE = ITEST                                                30.21
            IF (IOUTES .GT. ITEST) ITEST = IOUTES

!           --- synchronize nodes

            CALL SWSYNC                                                         40.30

!           --- carry out the output requests

!TIMG            CALL SWTSTA(9)                                                40.30
            CALL SWOUTP ( AC2   , SPCSIG, SPCDIR, COMPDA, XYTST ,         40.31
     &                    KGRPNT, XCGRID, YCGRID, KGRBND, OURQT )         40.51 40.31
!TIMG            CALL SWTSTO(9)                                                40.30
!
            IF (ERRPTS.GT.0) REWIND(ERRPTS)                               30.50
            ITEST = SAVITE                                                30.21


!           --- update time

            IF (NSTATM.EQ.1) THEN                                         40.00
              IF (NSTATC.EQ.1 .AND. IT.LT.MTC) TIMCO = TIMCO + DT         40.00
              CHTIME = DTTIWR(ITMOPT, TIMCO)                              40.00
              IF (NSTATC.EQ.1) WRITE (PRINTF, 222) CHTIME, TIMCO          40.00
 222          FORMAT(' Time of computation ->  ',A,' in sec:', F12.0)     40.00
            ENDIF                                                         40.00

         RETURN
         END   

