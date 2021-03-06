C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_calc.F,v 1.3 2005/04/06 18:45:20 jmc Exp $
C $Name:  $

#include "PP81_OPTIONS.h"

CBOP
C !ROUTINE: PP81_CALC

C !INTERFACE: =======================================================
       subroutine PP81_CALC(
     I     bi, bj, myTime, myThid )

C !DESCRIPTION: \bv
C     /==========================================================\
C     | SUBROUTINE PP81_CALC                                     |
C     | o Compute all PP81 fields defined in PP81.h              |
C     |==========================================================|
C     | This subroutine is based on SPEM code                    |
C     \==========================================================/
      IMPLICIT NONE
C
C--------------------------------------------------------------------

C JMK
C global parameters updated by pp_calc
C     PPviscAz   - PP eddy viscosity coefficient              (m^2/s)
C     PPdiffKzT  - PP diffusion coefficient for temperature   (m^2/s)
C
C \ev

C !USES: ============================================================
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"
#include "DYNVARS.h"
#include "PP81.h"
#include "FFIELDS.h"
#include "GRID.h"
#ifdef ALLOW_AUTODIFF_TAMC
#include "tamc.h"
#include "tamc_keys.h"
#else /* ALLOW_AUTODIFF_TAMC */
      integer ikppkey
#endif /* ALLOW_AUTODIFF_TAMC */

C !INPUT PARAMETERS: ===================================================
c Routine arguments
c     bi, bj - array indices on which to apply calculations
c     myTime - Current time in simulation

       INTEGER bi, bj
       INTEGER myThid
       _RL     myTime

#ifdef ALLOW_PP81

C !LOCAL VARIABLES: ====================================================
c Local constants
C     imin, imax, jmin, jmax  - array computation indices
C     RiNumber                - Richardson Number
      INTEGER I, J, K, Km1, JJ
      INTEGER   iMin ,iMax ,jMin ,jMax,di
      _RL     denom, PPviscTmp, buoyFreq,rhot,tempu,tempv,oldK,Ri
      _RL     b0, buoyFreqf, buoyFreqc, PPviscold,zsum,zsums
      _RL     rhoS(1:Nr),RS(1:Nr)
      _RL     dzc,dzp,ec,ep,es,epss(-1:0),epsw(-1:0),dz
!      _RL     bF(1:Nr)
!      _RL     theta_mcb(1:Nr),theta_mcb3(1:Nr)
C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
      CHARACTER*(1024) msgBuf

CEOP
      iMin = 2-OLx
      iMax = sNx+OLx-1
      jMin = 2-OLy
      jMax = sNy+OLy-1

C      DO K = 2, Nr

C ---------------------------------------------------------------------------
C replace theta with theta_mcb
C make linear temp profile
C ---------------------------------------------------------------------------
!       DO K=1,Nr
!         theta_mcb(K) = 20./3000. * rC(K) + 30.
!       ENDDO
!       theta_mcb3 = theta_mcb

C reverse surface layers
!       DO K=1,10
!          theta_mcb3(K+10) = theta_mcb(K)
!          theta_mcb3(K)    = theta_mcb(K+10)
!       ENDDO
C only reverse first layer
!       theta_mcb3(1+1) = theta_mcb(1)
!       theta_mcb3(1)   = theta_mcb(1+1)
!       theta_mcb= theta_mcb3 !switch

C Code by Jody Klymak has changed
C possible bugs
C 1) rhot = -tAlpha*(theta(I,J,K,bi,bj))*rhoNil
C K is not assigned and rho0 = 1024.0 is missing (?) This gives a  very low surface density 
C => the surface is always stable => causing a really high buoyancy freq (buoyFreq) 
C => giving a really high viscosity K in the case of near surface overturns.
C It almost seems if this is all purposely done, but I do not know why?
C 2) then something goes wrong with the sorting loop (I did not nail it), 
C because DZ is not equal over the sorted overturn 
C => -90 m at the surface and 100 m in the lower layer => gives unequal K
C 3) this incorrect soring causes a jump in the vertical density gradient, 
C which should be constant => gives a jump in N => gives a jump in K C
C 
C K will be 0.2*dz^2 N where dz is the displacement of water from its
C sorted depth and N is the brunt-vaisala frequency of the sorted profile.

      DO J=jMin,jMax
         DO I=iMin,iMax
            K=1
            CALL FIND_RHO_SCALAR(theta(I,J,K,bi,bj), salt(I,J,K,bi,bj),
     $           totPhiHyd(I,J,K,bi,bj),rhot,myThid )
            rhoS(1)=rhot
            RS(1)=rC(1)
            
            PPeps(I-1,J-1,1,bi,bj)=0.0
C           eps(k-1) = (dz(k-1)*eps0(k-1) +dz(k)*eps0(k))/(dz(k-1)+dz(k)) 
            ep = 0.0
            dzp = 0.0

            PPviscAr(I,J,1,bi,bj) = viscArNr(1)
            PPviscold = PPviscAr(I,J,1,bi,bj) ! at previous cell center

            DO K=2,Nr
               CALL FIND_RHO_SCALAR(theta(I,J,K,bi,bj), salt(I,J,K,bi
     $              ,bj),totPhiHyd(I,J,K,bi,bj),rhot,myThid )
               rhoS(K)=rhot
               RS(K)=rC(K)
C$$$               WRITE(msgBuf, '(A,I10.10,A,E10.4,A,E10.4)') 'Hellok ', K
C$$$     $              -1,' ',theta(I,J,K,bi,bj),' ',rhot
C$$$               CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
C$$$     &              SQUEEZE_RIGHT , 1)

               IF ( (rhoS(K).LT.rhoS(K-1)).AND.(maskC(I,J,K,bi
     &              ,bj).GT.0)) THEN
                  JJ=K-1
                  DO WHILE ( (JJ.GT.0).AND.(rhoS(K).LT.rhoS(JJ)) )
C                     write(*,*) K,JJ,rhoS(K),rhoS(JJ)   
                     JJ=JJ-1
                  ENDDO
                  rhoS(JJ+1:K)=cshift(rhoS(JJ+1:K),-1)
                  RS(JJ+1:K)=cshift(RS(JJ+1:K),-1)
               ENDIF
            ENDDO


C RS-R is dz....
C recip_drC=inverse distanance between centers, 
C first is between surface and first center
C diffKrNrT(K) = viscArNr(K) = background value   

            PPdiffKr(I,J,K,bi,bj) = MAX(PPviscAr(I,J,1,bi,bj),
     &              diffKrNrT(1))
C N at surface = zero or uses gradient
            b0 = MAX(-gravity*mass2rUnit*
     &              (rhoS(1) - rhoS(2))*recip_drC(2),0. _d 0)
C            b0 = 0.
            DO di=-1,0
               epss(di)=0.0
               epsw(di)=0.0
            ENDDO

            DO K=1,Nr
               IF (K.LT.Nr) THEN
                  buoyFreqf = -gravity*mass2rUnit*
     &              (rhoS(K) - rhoS(K+1))*recip_drC(K+1)
               ELSE
C N zero OR not zero near bottom (at the end of array)
                  buoyFreqf = -gravity*mass2rUnit*
     &              (rhoS(K-1) - rhoS(K))*recip_drC(K)
C                  buoyFreqf = 0.
               ENDIF 
               buoyFreqf = MAX(buoyFreqf,0. _d 0) ! not < 0
               buoyFreqc = (buoyFreqf + b0)*0.5   ! mean at cell center

C viscosity at cell center at K
               PPviscTmp = MAX( viscArNr(K), 0.2*(RS(K)-rC(K))*
     &                        (RS(K)-rC(K))*sqrt(buoyFreqc)) 

               IF (K.GT.1) THEN
                  Km1=K-1

C viscosity at cell face above center at K
                  PPviscAr(I,J,K,bi,bj) = MAX(0.5*(PPviscTmp+PPviscold),
     &                 viscArNr(K)) 
C Put an upper limit on viscosity to prevent instability when 
C explicit viscosity is C used (e.g. for nonhydrostatic case) SAL
                  PPviscAr(I,J,K,bi,bj) = MIN( 300. _d 0, 
     &                                  PPviscAr(I,J,K,bi,bj) )

                  PPdiffKr(I,J,K,bi,bj) = MAX( PPviscAr(I,J,K,bi,bj),
     &                                      diffKrNrT(K) )

C Compute Epsilon for diagnostics:
C 
C need to caclulate Im1 and Jm1 epsilon unfortunately...  Here at
C i-1,j-1 we average the west nu(du/dz)^2 at i-1 and i, and the south
C nu(dv/dv)^2 at j-1 and j, and then add them
C
C dz is calculated from the face distances, with the cells assumed to be
C half way.  Note the use of hfacW and hfacS to make these correct near
C bathy.  
                  zsum=0.
                  ec=0.0
                  zsums=0.
                  es=0.
                  DO di=-1,0
                     IF (hfacW(I+di,J-1,K,bi,bj).GT.0.000001) THEN
                        dz = 0.5*(drF(K)*hfacW(I+di,J-1,K,bi,bj)
     $                       +drF(Km1)*hfacW(I+di,J-1,Km1,bi,bj))
                        IF (dz.GT.0.00001) THEN
                           tempu = (uVel(I+di,J-1,Km1,bi,bj)-uVel(I+di,J
     $                          -1,K,bi,bj))/dz
                           epsw(di)=tempu*tempu*PPviscAr(I+di,J-1,K,bi
     $                          ,bj)
                           ec=ec+epsw(di)*dz
                           zsum = zsum+dz
                        ENDIF
                     ELSE 
C                       This face is on the seafloor.  set epsilon=the
C                       previous and dz = half the face.
                        dz=0.5*(drF(Km1)*hfacW(I+di,J-1,Km1,bi ,bj))
                        ec=ec+epsw(di)*dz
                        zsum = zsum+dz
                     ENDIF
C Now do the v-component                  
                     IF (hfacS(I-1,J+di,K,bi,bj).GT.0.000001) THEN
                        dz = 0.5*(drF(K)*hfacS(I-1,J+di,K,bi,bj)
     $                       +drF(Km1)*hfacS(I-1,J+di,Km1,bi,bj))
                        IF (dz.GT.0.00001) THEN
                           tempu = (vVel(I-1,J+di,Km1,bi,bj)-vVel(I-1,J
     $                          +di,K,bi,bj))/dz
                           epss(di)=tempu*tempu*PPviscAr(I-1,J+di,K,bi
     $                          ,bj)
                           es = es+epss(di)*dz
                           zsums = zsums+dz
                        ENDIF
                     ELSE
C                       This face is on the seafloor.  set epsilon=the
C                       previous and dz = half the face.
                        dz=+0.5*(drF(Km1)*hfacS(I-1,J+di,Km1 ,bi,bj))
                        es = es+epss(di)*dz
                        zsums = zsums+dz
                     ENDIF
                  ENDDO
C                 take the average of the du/dz terms
                  IF (zsum.GT.0.00001) THEN
                     ec=ec/zsum
                  ENDIF
C                 take the average of the dv/dz terms
                  IF (zsums.GT.0.00001) THEN
                     es=es/zsums
                  ENDIF
C add the u and v dissipations:
                  ec=es+ec

C Note this ec is defined on cell faces K=2..NR at the center of the
C cells (i.e. at XC), so its above the density variables.
C
C So to get at the center of the cells, just average this one and the previous one.  And its a true average because the 

                  PPeps(I-1,J-1,Km1,bi,bj) = 0.5*(ep+ec)
                  IF (Km1.EQ.1) THEN
                     PPeps(I-1,J-1,Km1,bi,bj) = ec
                  ENDIF
                  ep=ec
               ENDIF
c$$$               WRITE(msgBuf, '(A,I10.10,A,E10.4,A,E10.4)') 'Hellok ', K
c$$$     $              -1,' ',theta(I,J,K,bi,bj),' ',PPeps(I-1,J-1,Km1,bi
c$$$     $              ,bj)
c$$$               CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
c$$$     &              SQUEEZE_RIGHT , 1)


               b0        = buoyFreqf ! at previous cell face
               PPviscold = PPviscTmp ! at previous cell center
            ENDDO
C           ENDDO K
C     set on K=Nr
            PPeps(I-1,J-1,Nr,bi,bj) =ep

         ENDDO
C           ENDDO J
      ENDDO
C           ENDDO I

#endif /* ALLOW_PP81 */

      RETURN
      END

