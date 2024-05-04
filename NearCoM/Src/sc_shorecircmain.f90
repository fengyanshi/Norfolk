SUBROUTINE SHORECIRC_INITIALIZATION
     USE GLOBAL
     IMPLICIT NONE
     INTEGER::ISTAGE

     CALL READ_INPUT

     CALL INDEX

! allocate variables
     CALL ALLOCATE_VARIABLES

     CALL INITIALIZATION


END SUBROUTINE SHORECIRC_INITIALIZATION



SUBROUTINE SHORECIRC_CYCLE
     USE GLOBAL



     IMPLICIT NONE
     INTEGER::ISTAGE

! update three variables
     Eta0=Eta
     Ubar0=Ubar
     Vbar0=Vbar  

     CALL UPDATE_MASK

     CALL EXCHANGE

! updates in exchange  

         CALL OneWayCoupling







     CALL ESTIMATE_DT(Mloc,Nloc,JS1,JS2,U,V,H,MinDepthFrc,DT,CFL,TIME)






!  Runge-Kutta Scheme 
     ! 3-ORDER RUNGE-KUTTA TIME STEPPING
     DO ISTAGE=1,3
 
       CALL FLUXES

       CALL BOUNDARY_CONDITION

       IF(WindForce)THEN
         CALL WIND_SPEED
       ENDIF

       CALL SOLVE_3D

       CALL SourceTerms   ! put sourceterms after fluxes in order to get eta_t

       CALL ESTIMATE_HUV(ISTAGE)  

! call boundary condition test
       CALL TIDE

       CALL EXCHANGE

! updates in exchange  

         CALL OneWayCoupling


       IF(SPONGE_ON)THEN
         CALL SPONGE_DAMPING
       ENDIF

     ENDDO
 
!   end Runge-Kutta Scheme





















! ANNA start




! end momentum balance

! ANNA end

END SUBROUTINE SHORECIRC_CYCLE
             

