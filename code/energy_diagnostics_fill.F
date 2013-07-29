C $Header: /u/gcmpack/MITgcm/pkg/diagnostics/diagnostics_fill_state.F,v 1.30 2006/02/16 23:03:16 molod Exp $
C $Name:  $

#include "DIAG_OPTIONS.h"

CBOP
C     !ROUTINE: DIAGNOSTICS_FILL_STATE
C     !INTERFACE:
      SUBROUTINE ENERGY_DIAGNOSTICS_FILL( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DIAGNOSTICS_FILL_STATE
C     | o Fill-in main code, state-variables diagnostics
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "NH_VARS.h"
#include "SURFACE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid     :: my Thread Id number
      INTEGER myThid 

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     == Local variables ==
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
      _RL uPbt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vPbt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uPbc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vPbc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uEbc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vEbc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uEbt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vEbt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Ebt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Ebc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL Conv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL N2(Nr)
      _RL uC,uS,uW,vC,vS,vW,wC,vN,vE
      _RL rhoC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoS,rhoT
      _RL U0E
      _RL PS,PW,W,EkS,EkW
      _RL U0C,U0S,U0N,V0C,V0W,V0E,V0N
      _RL ddz,kEpC
      _RL hkEpW,hkEpS,kEpS,kEpW
      _RL Ep,EpW,EpS,ptop,pbot,dqdz
      _RL upC,upE,upN,upS,upW,vpC,vpE,vpN,vpS,vpW
      _RL P0S(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL P0W(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL V0S(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL U0W(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL W1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL W2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL upupW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL upvpS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL upvpW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL vpvpS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ZS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL ZW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uE,uN
      INTEGER i,j,K,bi,bj,jj,ii
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      
C  Divide energy flux into various depth-integrated terms:
C  Linear baroclinic energy flux: u*P_BC
C  Linear barotropic energy flux: <u>*P_BT

      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO K=1,Nr
               DO j = 0,sNy+1    
                  DO i = 0,sNx+1
                     IF (K.EQ.1) THEN
                        P0W(i,j,bi,bj)=0.
                        P0S(i,j,bi,bj)=0.
                        U0W(i,j,bi,bj)=0. 
                        V0S(i,j,bi,bj)=0.
                        ZW(i,j,bi,bj)=0.000001
                        ZS(i,j,bi,bj)=0.000001
                     ENDIF
C     Get N2:  just do for the first i,j because its always the same.  N2 is defined on 
C     faces
                     IF ((j.EQ.1).AND.(i.EQ.1)) THEN
                        CALL FIND_RHO_SCALAR(
     I                       tRef(K), sRef(K), rC(K),
     O                       rhoT, myThid )
                        IF (K.EQ.1) THEN
                           CALL FIND_RHO_SCALAR(
     I                          tRef(K+1), sRef(K+1), rC(K+1),
     O                          rhoS, myThid )
                           N2(K)=gravity/rhoNil*(rhoS-rhoT)/drC(K)
                        ELSE 
                           CALL FIND_RHO_SCALAR(
     I                          tRef(K-1), sRef(K-1), rC(K-1),
     O                          rhoS, myThid )
                           N2(K)=gravity/rhoNil*(rhoT-rhoS)/drC(K)
                        ENDIF 
                     ENDIF
                     PS = 0.5*(totPhiHyd(i,j,K,bi,bj)+totPhiHyd(i,j-1,K
     $                    ,bi,bj)) 
                     PW = 0.5*(totPhiHyd(i-1,j,K,bi,bj)+totPhiHyd(i,j,K
     $                    ,bi,bj))
                    IF (nonHydroStatic) THEN
                      PS = PS+0.5*(phi_nh(i,j,K,bi,bj)+phi_nh(i,j-1,K
     $                   ,bi,bj)) 
                      PW = PS+0.5*(phi_nh(i-1,j,K,bi,bj)+phi_nh(i,j,K
     $                     ,bi,bj))
                    ENDIF
C     Depth mean pressure
                     P0W(i,j,bi,bj)=P0W(i,j,bi,bj)+PW*DRf(K)*
     $                    hFacW(i,j,K,bi,bj)
                     P0S(i,j,bi,bj)=P0S(i,j,bi,bj)+PS*DRf(K)*
     $                    hFacS(i,j,K,bi,bj)
C     Depth mean velocities
                     U0W(i,j,bi,bj)=U0W(i,j,bi,bj)+uVel(i,j,K,bi,bj)
     $                    *DRf(K)*hFacW(i,j,K,bi,bj)
                     V0S(i,j,bi,bj)=V0S(i,j,bi,bj)+vVel(i,j,K,bi,bj)
     $                    *DRf(K)*hFacS(i,j,K,bi,bj)
                     ZW(i,j,bi,bj)=ZW(i,j,bi,bj)+DRf(K)*
     $                    hFacW(i,j,K,bi,bj)
                     ZS(i,j,bi,bj)=ZS(i,j,bi,bj)+DRf(K)*
     $                    hFacS(i,j,K,bi,bj)
                     IF (K.EQ.NR) THEN
C Get the means!
                        P0W(i,j,bi,bj)=P0W(i,j,bi,bj)/ZW(i,j,bi,bj)
                        P0S(i,j,bi,bj)=P0S(i,j,bi,bj)/ZS(i,j,bi,bj)
                        U0W(i,j,bi,bj)=U0W(i,j,bi,bj)/ZW(i,j,bi,bj)
                        V0S(i,j,bi,bj)=V0S(i,j,bi,bj)/ZS(i,j,bi,bj)
                     ENDIF
C               Enddo: i,j,k
                  ENDDO
               ENDDO
            ENDDO
C done barotropic loop.  Now do baroclinic loop:

            DO K=1,Nr
               CALL FIND_RHO_2D(
     I              1, sNx, 1, sNy, K,
     I              theta(1-OLx,1-OLy,K,bi,bj), 
     I              salt(1-OLx,1-OLy,K,bi,bj),
     O              rhoC,
     I              K,bi,bj,myThid )
C               CALL FIND_RHO_SCALAR(theta(i,j,K,bi,bj), 
C     &              salt(i,j,K,bi,bj),
C     &              totPhiHyd(i,j,K,bi,bj),rhoC,myThid )

               DO j = 0,sNy+1    
                  DO i = 0,sNx+1
C W1 and W2 needed for next loop.  Eq 5.42
                     IF (K.eq.1) THEN
                        W1(i,j,bi,bj) = (ZW(i+1,j,bi,bj)*U0W(i+1,j,bi
     $                       ,bj)-ZW(i,j,bi,bj)*U0W(i,j,bi ,bj))
     $                       *recip_dxF(i,j,bi,bj)
C                        W1(i,j,bi,bj) = W1(i,j,bi,bj)+ (ZS(i,j+1,bi,bj)
C     $                       *V0S(i,j+1,bi,bj)- ZS(i,j,bi,bj)*V0S(i,j,bi
C     $                       ,bj))*recip_dyF(i,j,bi,bj)
                        W2(i,j,bi,bj) = (U0W(i+1,j,bi
     $                       ,bj)- U0W(i,j,bi,bj))*recip_dxF(i,j,bi,bj)
C                        W2(i,j,bi,bj) = W2(i,j,bi,bj)+ (V0S(i,j+1,bi,bj)
C     $                       - V0S(i,j,bi,bj))*recip_dyF(i,j,bi,bj)
                     ENDIF
C just in case:
                     rhoC(i,j)=rhoC(i,j)-rhoNil+rhoConst
c BT W Eq 5.42 K11
                     W = -W1(i,j,bi,bj)-W2(i,j,bi,bj)*rC(K)
C get fields at center, south face and west face of each cell, east face
C and north face.  For uS, vW etc these need to be calculated as
C averages of the appropriate surrounding cells.  In order to take into
C account cells that are partially full or even empty, we weight by
C hfacW for u and hfacS for v.
  
                     uW = uVel(i,j,K,bi,bj)
                     vS = vVel(i,j,K,bi,bj)
                     uE = uVel(i+1,j,K,bi,bj)
                     vN = vVel(i,j+1,K,bi,bj)
C  uC                    
                     ddz = 0.
                     uC=0.
                     uC = uW*hfacW(i,j,K,bi,bj)
                     ddz = ddz +hfacW(i,j,K,bi,bj) 
                     uC = uC + uE*hfacW(i+1,j,K,bi,bj)
                     ddz = ddz +hfacW(i+1,j,K,bi,bj) 
                     if (ddz.GT.0.) THEN
                        uC = uC/ddz
                     ENDIF
C vC
                     ddz = 0.
                     vC=0.
                     vC = vS*hfacS(i,j,K,bi,bj)
                     ddz = ddz +hFacS(i,j,K,bi,bj) 
                     vC = vC + vN*hfacS(i,j+1,K,bi,bj)
                     ddz = ddz +hfacS(i,j+1,K,bi,bj) 
                     if (ddz.GT.0.) THEN
                        vC = vC/ddz
                     ENDIF
C  uN 
                     ddz = 0.
                     uN=0.
                     DO jj = 0,1
                        DO ii = 0,1
                           uN = uN+uVel(i+ii,j+jj,K,bi,bj)*hfacW(i+ii,j
     $                          +jj,K,bi,bj)
                           ddz = ddz +hfacW(i+ii,j+jj,K,bi,bj) 
                        ENDDO
                     ENDDO
                     if (ddz.GT.0.) THEN
                        uN = uN/ddz
                     ENDIF
C  uS 
                     ddz = 0.
                     uS=0.
                     DO jj = -1,0
                       DO ii = 0,1
                          uS = uS+uVel(i+ii,j+jj,K,bi,bj)*hfacW(i+ii,j
     $                         +jj,K,bi,bj)
                          ddz = ddz +hfacW(i+ii,j+jj,K,bi,bj) 
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      uS = uS/ddz
                    ENDIF
C  vW 
                    ddz = 0.
                    vW=0.
                    DO jj = 0,1
                       DO ii = -1,0
                          vW = vW+vVel(i+ii,j+jj,K,bi,bj)*hfacS(i+ii,j
     $                         +jj,K,bi,bj)
                          ddz = ddz +hfacS(i+ii,j+jj,K,bi,bj) 
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      vW = vW/ddz
                    ENDIF
C  vE
                    ddz = 0.
                    vE=0.
                    DO jj = 0,1
                       DO ii = 0,1
                          vE = vE+vVel(i+ii,j+jj,K,bi,bj)*hfacS(i+ii,j
     $                         +jj,K,bi,bj)
                          ddz = ddz +hfacS(i+ii,j+jj,K,bi,bj) 
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      vE = vE/ddz
                    ENDIF

C depth-weighted mean of barotropic quantities.  
                    IF (rlowW(i,j,bi,bj).LT.0.) THEN
                       U0C = (-rlowW(i,j,bi,bj)*U0W(i,j,bi,bj)-rlowW(i+1
     $                      ,j,bi,bj)*U0W(i+1,j,bi,bj))/(-rlowW(i,j,bi
     $                      ,bj)-rlowW(i+1,j,bi,bj))
                    ELSE
                       U0C=0.
                    ENDIF
                    U0E = U0W(i+1,j,bi,bj)
C U0S:
                    U0S = 0.0
                    ddz = 0.
                    DO jj = -1,0
                       DO ii = 0,1
                          U0S = U0S+U0W(i+ii,j+jj,bi,bj)*(-rlowW(i+ii,j
     $                         +jj,bi,bj))
                          ddz = ddz -rlowW(i+ii,j+jj,bi,bj)
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      U0S = U0S/ddz
                    ENDIF
C U0N:
                    U0N = 0.0
                    ddz = 0.
                    DO jj = 0,1
                       DO ii = 0,1
                          U0N = U0N+U0W(i+ii,j+jj,bi,bj)*(-rlowW(i+ii,j
     $                         +jj,bi,bj))
                          ddz = ddz-rlowW(i+ii,j+jj,bi,bj)
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      U0N = U0N/ddz
                    ENDIF
C V0C
                    V0C = 0.0
                    ddz = 0.
                    DO jj = 0,1
C                       DO ii = 0,0
                       V0C = V0C+V0S(i,j+jj,bi,bj)*(-rlowS(i,j+jj
     $                      ,bi,bj))
                          ddz = ddz-rlowS(i,j+jj,bi,bj)
C                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      V0C = V0C/ddz
                    ENDIF
C V0W
                    V0W = 0.0
                    ddz = 0.
                    DO jj = 0,1
                       DO ii = -1,0
                          V0W = V0W+V0S(i+ii,j+jj,bi,bj)*(-rlowS(i+ii,j
     $                         +jj,bi,bj))
                          ddz = ddz-rlowS(i+ii,j+jj,bi,bj)
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      V0W = V0W/ddz
                    ENDIF
C V0E
                    V0E = 0.0
                    ddz = 0.
                    DO jj = 0,1
                       DO ii = 0,1
                          V0E = V0E+V0S(i+ii,j+jj,bi,bj)*(-rlowS(i+ii,j
     $                         +jj,bi,bj))
                          ddz = ddz -rlowS(i+ii,j+jj,bi,bj)
                       ENDDO
                    ENDDO
                    if (ddz.GT.0.) THEN
                      V0E = V0E/ddz
                    ENDIF
C V0N
                    V0N=V0S(i,j+1,bi,bj)
C primed quantities:
                    upW = uW-U0W(i,j,bi,bj)
                    upC = uC-U0C
                    upE = uE-U0E
                    upN = uN-U0N
                    upS = uS-U0S
                    vpS = vS-V0S(i,j,bi,bj)
                    vpN = vN-V0N
                    vpE = vE-V0E
                    vpW = vW-V0W
                    vpC = vC-V0C
C w in the center
                    IF (K.LT.Nr) THEN
                       wC = 0.5*(wVel(i,j,K,bi,bj)+wVel(i,j,K+1,bi,bj))
C Lets not be so precise w. W and then we need not worry about masks
C                       wS = 0.25*(wVel(i,j,K,bi,bj)+wVel(i,j,K+1,bi,bj)
C    $                    +wVel(i,j-1,K,bi,bj)+wVel(i,j-1,K+1,bi,bj))
C                       wW = 0.25*(wVel(i,j,K,bi,bj)+wVel(i,j,K+1,bi,bj)
C    $                    +wVel(i-1,j,K,bi,bj)+wVel(i-1,j,K+1,bi,bj))
                    ELSE
                       wC = wVel(i,j,K,bi,bj)
                    ENDIF
C                   ---  zero out the arrays
                    if (K.EQ.1) THEN
                       uPbc(i,j,bi,bj)=0.
                       vPbc(i,j,bi,bj)=0.
                       uEbc(i,j,bi,bj)=0.
                       vEbc(i,j,bi,bj)=0.
                       Ebc(i,j,bi,bj)=0.
                       Conv(i,j,bi,bj)=0.
                       uPbt(i,j,bi,bj)=0.
                       vPbt(i,j,bi,bj)=0.
                       upupW(i,j,bi,bj)=0.
                       upvpW(i,j,bi,bj)=0.
                       upvpS(i,j,bi,bj)=0.
                       vpvpS(i,j,bi,bj)=0.
                    ENDIF
C Pressure U*P. 
                    PS = 0.5*(totPhiHyd(i,j,K,bi,bj)+totPhiHyd(i,j-1,K
     $                   ,bi,bj)) 
                    PW = 0.5*(totPhiHyd(i-1,j,K,bi,bj)+totPhiHyd(i,j,K
     $                   ,bi,bj))
                    IF (nonHydroStatic) THEN
                      PS = PS+0.5*(phi_nh(i,j,K,bi,bj)+phi_nh(i,j-1,K
     $                   ,bi,bj)) 
                      PW = PS+0.5*(phi_nh(i-1,j,K,bi,bj)+phi_nh(i,j,K
     $                     ,bi,bj))
                    ENDIF
C               Baroclinic flux:
                    uPbc(i,j,bi,bj)=uPbc(i,j,bi,bj)+(PW-P0W(i,j,bi,bj))
     $                   *(upW)*DRf(K)*hFacW(i,j,K,bi,bj)
                    vPbc(i,j,bi,bj)=vPbc(i,j,bi,bj)+(PS-P0S(i,j,bi,bj))
     $                   *(vpS)*DRf(K)*hFacS(i,j,K,bi,bj)
C               Baroclinic KE:
                    kEpC = 0.5*(upC*upC+vpC*vpC + wC*wC)
                    kEpW = 0.5*(upW*upW+vpW*vpW + wC*wC)
                    kEpS = 0.5*(upS*upS+vpS*vpS + wC*wC)
C eqn 5.51:
                    hkEpW = U0W(i,j,bi,bj)*upW+V0W*vpW
                    hkEpS = U0S*upS+V0S(i,j,bi,bj)*vpS
C Potential energy!!! Here we will use the linear approximation under the assumption that
C the stratification is not too non-WKB
                    Ep = gravity*gravity*rhoC(i,j)*rhoC(i,j) 
                    Ep=Ep/2/rhoNil/rhoNil/N2(K)
C Again a shortcut since these only appear in the advective terms...
                    EpW = Ep
                    EpS = Ep
C               Advective flux:
                    uEbc(i,j,bi,bj) = uEbc(i,j,bi,bj)+uW*(kEpW+hkEpW
     $                   +EpW)*DRf(K)*hFacW(i,j,K,bi,bj)
                    vEbc(i,j,bi,bj) = vEbc(i,j,bi,bj)+vS*(kEpS+hkEpS
     $                   +EpS)*DRf(K)*hFacS(i,j,K,bi,bj)
C KEp and APEp sums:
                    Ebc(i,j,bi,bj)=Ebc(i,j,bi,bj)+kEpC*DRf(K)*hFacC(i,j
     $                   ,K,bi,bj)
                    Ebc(i,j,bi,bj)=Ebc(i,j,bi,bj)+Ep*DRf(K)*hFacC(i,j,K
     $                   ,bi,bj)
C Conversion (5.102)
                    Conv(i,j,bi,bj)=Conv(i,j,bi,bj)+rhoC(i,j)*gravity*W
     $                   *DRf(K)*hFacC(i,j,K,bi,bj)/rhoNil

                    IF ((i.eq.100).AND.(j.eq.1)) THEN
                       WRITE(msgBuf,
     $                      '(A,I2,A,E10.3,A,E10.4,A,E10.4,A,E10.4)')
     $                      'Hellok ', K,'rho: ',rhoC(i,j),' Ep: ',Ep,'
     $                      N2: ' ,N2(K),' Ke: ',kEpC
C                       CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
C     $                      SQUEEZE_RIGHT , 1)
                       WRITE(msgBuf,
     $                      '(A,I2,A,E10.3,A,E10.4,A,E10.4,A,E10.4)')
     $                      'Helloj ',j,' vC: ',vC,' hkEpW: ',hkEpW,'
     $                      EpW: ' ,EpW,' kEpW: ',kEpW
C                       CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
C     $                      SQUEEZE_RIGHT , 1)
                       WRITE(msgBuf,
     $                      '(A,I2,A,E10.3,A,E10.4,A,E10.4,A,E10.4)')
     $                      'Helloi ',i,' W: ',W ,' Conv: ' ,Conv(i,j,bi
     $                      ,bj),' W1: ' ,W1(i,j,bi,bj) ,' W2: ',W2(i,j
     $                      ,bi,bj)
C                       CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
C     $                      SQUEEZE_RIGHT , 1)
                    ENDIF

C Need a dq/dz term.  Ughhhh....  and phi_nh is at the grid centers, so need a 
C centered difference.  
                    IF (nonHydroStatic) THEN
                      IF (K.EQ.1) THEN
                        ptop=phi_nh(i,j,K,bi,bj)
                        pbot = phi_nh(i,j,K+1,bi,bj)
                        dqdz = (ptop-pbot)/dRC(1)
                      ELSE IF ((K.EQ.NR)) THEN
                        ptop=phi_nh(i,j,K-1,bi,bj)
                        pbot = phi_nh(i,j,K,bi,bj)
                        dqdz = (ptop-pbot)/dRC(NR)
                      ELSE
                         ptop = hFacC(i,j,K-1,bi,bj)*phi_nh(i,j,K-1,bi
     $                        ,bj) 
                         ptop = ptop + hFacC(i,j,K,bi,bj)*phi_nh(i,j,K
     $                        ,bi,bj) 
                         ddz = hFacC(i,j,K-1,bi,bj)+hFacC(i,j,K,bi,bj)
                         IF (ddz>0.) THEN
                            ptop = ptop/ddz
                         ELSE 
                            ptop = 0.0
                         ENDIF
                         pbot = hFacC(i,j,K+1,bi,bj)*phi_nh(i,j,K+1,bi
     $                        ,bj) 
                         pbot = pbot + hFacC(i,j,K,bi,bj)*phi_nh(i,j,K
     $                        ,bi,bj) 
                         ddz = hFacC(i,j,K,bi,bj)+hFacC(i,j,K+1,bi,bj)
                         if (ddz>0.) THEN
                            pbot = pbot/ddz
                         ELSE
                            pbot=0.0
                         ENDIF
                         dqdz = (ptop-pbot)/DRf(NR)
                      ENDIF
                      Conv(i,j,bi,bj)=Conv(i,j,bi,bj)-dqdz*W *DRf(K)
     $                     *hFacC(i,j,K,bi,bj)
                    ENDIF
C Aho: This is a barotropic term, but we need to average them now.
C  5.79, 5.80, 5.86  
                    upupW(i,j,bi,bj)= upupW(i,j,bi,bj)+upW*upW*DRf(k)
     $                   *hFacW(i,j,K,bi,bj)
                    upvpS(i,j,bi,bj)=upvpS(i,j,bi,bj)+upS*vpS*DRf(k)
     $                   *hFacS(i,j,K,bi,bj)
                    upvpW(i,j,bi,bj)=upvpW(i,j,bi,bj)+upW*vpW*DRf(k)
     $                   *hFacW(i,j,K,bi,bj)
                    vpvpS(i,j,bi,bj)=vpvpS(i,j,bi,bj)+vpS*vpS*DRf(k)
     $                   *hFacS(i,j,K,bi,bj)
C This is the barotropic:
                    uPbt(i,j,bi,bj)=uPbt(i,j,bi,bj)+P0W(i,j,bi,bj)*uW
     $                   *DRf(K)*hFacW(i,j,K,bi,bj)
C There should be a mean q (non hydrostatic) term as well.  
                    vPbt(i,j,bi,bj)=vPbt(i,j,bi,bj)+P0S(i,j,bi,bj)*vS
     $                   *DRf(K)*hFacS(i,j,K,bi,bj)
                 ENDDO 
              ENDDO
           ENDDO
C Need to add to the conversion term the Ah0 term ( 5.86, 5.79, and 5.80)
           DO j = 1,sNy    
              DO i = 1,sNx
                 U0C=0.
                 IF (rlowW(i,j,bi,bj).LT.0.) THEN
                    U0C = (-rlowW(i,j,bi,bj)*U0W(i,j,bi,bj)-rlowW(i+1,j
     $                   ,bi,bj)*U0W(i+1,j,bi,bj))/(-rlowW(i,j,bi,bj)
     $                   -rlowW(i+1,j,bi,bj))
                 ENDIF
                 V0C = 0.0
                 ddz = 0.
                 DO jj = 0,1
                    ii=0
C                    DO ii = 0
                    V0C = V0C+V0S(i+ii,j+jj,bi,bj)*(-rlowS(i+ii,j+jj ,bi
     $                   ,bj))
                       ddz = ddz -rlowS(i+ii,j+jj,bi,bj)
C                    ENDDO
                 ENDDO
                 if (ddz.GT.0.) THEN
                    V0C = V0C/ddz
                 ENDIF
                 Conv(i,j,bi,bj)=Conv(i,j,bi,bj)+U0C*(upupW(i+1,j,bi
     $                ,bj)-upupW(i,j,bi,bj))/dxF(i,j,bi,bj)
                 Conv(i,j,bi,bj)=Conv(i,j,bi,bj)+U0C*(upvpS(i,j+1,bi
     $                ,bj)-upvpS(i,j,bi,bj))/dyF(i,j,bi,bj)
                 Conv(i,j,bi,bj)=Conv(i,j,bi,bj)+V0C*(upvpW(i+1,j,bi
     $                ,bj)-upvpW(i,j,bi,bj))/dxF(i,j,bi,bj)
                 Conv(i,j,bi,bj)=Conv(i,j,bi,bj)+V0C*(vpvpS(i,j+1,bi
     $                ,bj)-vpvpS(i,j,bi,bj))/dyF(i,j,bi,bj)
C Do the barotroic terms here too...
                 Ebt(i,j,bi,bj) = -0.5*(U0C*U0C+V0C*V0C)*R_low(i,j,bi
     $                ,bj)
                 Ebt(i,j,bi,bj)= Ebt(i,j,bi,bj)+0.5*etaH(i,j,bi,bj)
     $                *etaH(i,j,bi,bj)*gravity
                 EkS = -0.5*(U0C*U0C+V0S(i,j,bi,bj)*V0S(i,j,bi,bj))
     $                *rlowS(i,j,bi,bj)
                 EkW = -0.5*(U0W(i,j,bi,bj)*U0W(i,j,bi,bj)+V0C*V0C)
     $                *rlowW(i,j,bi,bj)
                 uEbt(i,j,bi,bj) = U0W(i,j,bi,bj)*EkW
                 vEbt(i,j,bi,bj) = V0S(i,j,bi,bj)*EkS
                 IF ((i.eq.100).AND.(j.eq.1)) THEN
C                    WRITE(msgBuf,
C     $                   '(A,I2,A,E10.3,A,E10.4,A,E10.4,A,E10.4)') 'i ',
C     $                   i,'Ebt: ',Ebt(i,j,bi,bj),' U0W: ',U0W(i,j,bi
C     $                   ,bj),' R_low: ' ,R_low(i,j,bi,bj),' EkW: ',EkW
C                    CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
C     $                   SQUEEZE_RIGHT , 1)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
      ENDDO
C BT
      CALL DIAGNOSTICS_FILL(Ebt, 'SDIAG1  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(uPbt,'SDIAG2  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(vPbt,'SDIAG3  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(uEbt,'SDIAG4  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(vEbt,'SDIAG5  ',0,1,0,1,1,myThid)
C BC
      CALL DIAGNOSTICS_FILL(Ebc, 'SDIAG6  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(uPbc,'SDIAG7  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(vPbc,'SDIAG8  ',0,1,0,1,1,myThid)
      CALL DIAGNOSTICS_FILL(uEbc,'SDIAG9  ',0,1,0,1,1,myThid)
C      CALL DIAGNOSTICS_FILL(vEbc,'SDIAG10 ',0,1,0,1,1,myThid)
C Conversion
      CALL DIAGNOSTICS_FILL(Conv,'SDIAG10 ',0,1,0,1,1,myThid)
C set up in data.diagnostics:
C 
       
             

#endif /* ALLOW_DIAGNOSTICS */
      
      RETURN 
      END
