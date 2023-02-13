*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST
      PARAMETER (NPOINT = 77)
      DIMENSION EN(NPOINT), SPEC(NPOINT)
*
* hall_d_gamma_spectrum_on_target.vect file from Pavel
*
      SAVE LFIRST, ENMIN, ENMAX
      DATA LFIRST / .TRUE. /

      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*  |  *** User initialization ***
         LUNRD = NINT(WHASOU(1))
         WRITE(LUNOUT,*)
         WRITE(LUNOUT,'(A,132A)') ("*",I=1,132)
         WRITE(LUNOUT,*)
         WRITE(LUNOUT,*)
         WRITE(LUNOUT,'(A,A)') "     Solid Spectrum",  
     &   " ( file spectrum_neutron_coil.dat )"
         ENMIN = 1.0D-5
         ENMAX = 1.05D0
         WRITE(LUNOUT,*)
         WRITE(LUNOUT,*)
         WRITE(LUNOUT,'(A,132A)') ("*",I=1,132)
         WRITE(LUNOUT,*)
*        Read the energies and corresponding spectrum values
         DO 1 I = 1, NPOINT    
            READ(LUNRD,*,END=999) E1, E2, SPEC(I)
            EN(I) = (E1+E2)*HLFHLF
            IF (I .EQ. 1) THEN
               EN(I) = E1
            END IF
            IF (I .EQ. NPOINT) THEN
               EN(I) = E2
            END IF
*        Energy is written in GeV in the file, and here I needed in GeV
 1       CONTINUE
 999     CONTINUE
*        First and last points extrapolated
         ENMIN = EN(1)
         ENMAX = EN(NPOINT)
*         SPEC(1) = 0.700503
  
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
      RAN01 = FLRNDM(DUMMY)
*     Sample an energy from a flat distribution between Enmin and Enmax
      ENERGY = ENMIN + (ENMAX - ENMIN) * RAN01
*     Interpolate the spectrum value and set the weight equal to it
      DO 2 I = 1, NPOINT - 1
         IF(ENERGY .GT. EN(I) .AND. ENERGY .LE. EN(I+1)) THEN
            WEIGHT = SPEC(I) + (SPEC(I+1) - SPEC(I)) *
     &                         (ENERGY - EN(I))/(EN(I+1) - EN(I))
            GO TO 3
         END IF
 2    CONTINUE
      STOP "SOURCE routine: Energy interval not found"
 3    CONTINUE
* Npflka is the stack counter: of course any time source is called it
* must be =0
      NPFLKA = NPFLKA + 1
* Wt is the weight of the particle
      WTFLK  (NPFLKA) = WEIGHT
      WEIPRI = WEIPRI + WTFLK (NPFLKA)
* Particle type (1=proton.....). Ijbeam is the type set by the BEAM
* card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
      END IF
*  |
*  +-------------------------------------------------------------------*
* From this point .....
* Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
* User dependent flag:
      LOUSE  (NPFLKA) = 0
* User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
* User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE
* Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
* ... to this point: don't change anything
* Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO
* Group number for "low" energy neutrons, set to 0 anyway
      IGROUP (NPFLKA) = 0
* Kinetic energy of the particle (GeV)
      TKEFLK (NPFLKA) = ENERGY
* Particle momentum
*     For photons, momentum = energy
      PMOFLK (NPFLKA) = ENERGY
* Cosines (tx,ty,tz)
      TXFLK  (NPFLKA) = UBEAM
      TYFLK  (NPFLKA) = VBEAM
      TZFLK  (NPFLKA) = WBEAM
* Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER
*     Particle coordinates
      RAN01 = FLRNDM(DUMMY)
      XFLK   (NPFLKA) = XBEAM -0.45+RAN01*0.9
      RAN01 = FLRNDM(DUMMY)
      ZFLK   (NPFLKA) = ZBEAM-50.0+RAN01*100
      YFLK   (NPFLKA) = YBEAM
*      WRITE(LUNOUT,*) 'X', XFLK (NPFLKA)
*      WRITE(LUNOUT,*) 'Y', YFLK (NPFLKA)
*      WRITE(LUNOUT,*) 'Z', ZFLK (NPFLKA)

*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
*  Flag this is prompt radiation
      LRADDC (NPFLKA) = .FALSE.
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
      END
