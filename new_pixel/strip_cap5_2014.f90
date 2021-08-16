!       ****************************************************** 
!            LAPLACE SOLVER FOR FINITE DEPTH 
!            STRUCTURE OF 5 PARALLEL STRIPS
!	     INCLUDING THE CONTRIBUTION OF AIR.
!            INPUT IS ON THE FILE instrip.2d AND
!              OUTPUT IS FILE strip.2d.out 
!                 2-D Calculation.
!                      May 1995
!       ******************************************************

!       ******************************************************
!                 DEFINITIONS OF VARIABLES
!       ******************************************************

        INTEGER DIM
        PARAMETER(DIM=262144)
        REAL DATA(2*DIM),P(DIM,2),W(DIM)
        REAL P1(2*DIM),PS(2*DIM),PSS(DIM,2),P1S(DIM,2)
        REAL R,PI,C00,C01,C02,EAIR,ESI,C1D
        REAL LM,D,C1,C2,C3,C4,C5,C0
        INTEGER NN,RUNS,Z,LJ
        INTEGER RES,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,PITCH,L,CYCLE
	COMMON DATA
        RES=DIM
	ESI=11.9
	EAIR=1.
        OPEN(UNIT=10,FILE='instrip.2d')
    	OPEN(UNIT=22,FILE='strip.2d.out')
        READ(10,*)RUNS
        DO 1968,Z=1,RUNS

!       *****************************************************
!               INPUT FROM FILE  instrip.2d
!       ***************************************************** 

        READ(10,*)RES,LJ,PITCH,D

!       ****************************************************
!               HYPOTHETICAL SOLUTION ON SURFACE
!       ****************************************************

        PI=3.14159
        CYCLE=0
   
	A5=(RES-LJ)/2
	A4=A5-PITCH+1
	A3=A4-LJ-1
	A2=A3-PITCH+1
	A1=A2-LJ-1
	A6=(RES+LJ)/2+1
	A7=A6+PITCH-1
	A8=A7+LJ+1
 	A9=A8+PITCH-1
 	A10=A9+LJ+1
        print *, 'placement of strips on the grid'
        print *, 'A1=', A1
        print *, 'A2=', A2
        print *, 'A3=', A3
        print *, 'A4=', A4
        print *, 'A5=', A5
        print *, 'A6=', A6
        print *, 'A7=', A7
        print *, 'A8=', A8
        print *, 'A9=', A9
        print *, 'A10=',A6
!    **************************************************
!       Initialize arrays:
!       Set the voltage on the central strip equal to 1.
!       and zeros everywhere else in the Voltage arrays and
!       the Fourier arrays
!    ********************************************************
        DO 1,K=1,RES
        P(K,2)=0.0
	    P(K,1)=0.
        IF((K.GT.A5).AND.(K.LT.A6))THEN
        P(K,1)=1.0
        END IF
  1     CONTINUE
  
        DO 3,K=A4-1,A5+1
        P(K,1)=(REAL(K)-REAL(A4-1))/(REAL(A5+1)-REAL(A4-1))
  3     CONTINUE
        DO 5,K=A6-1,A7+1
        P(K,1)=(REAL(A7+1)-REAL(K))/(REAL(A7+1)-REAL(A6-1))
  5     CONTINUE
        NN=RES

!       ****************************************************       
!          HYPOTHETICAL SOLUTION ON SURFACE NOW IS SET
!       ****************************************************
        
 222    CONTINUE        
        CYCLE=CYCLE+1

!       ***************************************************
!           FOURIER TRANSFORM OF SURFACE POTENTIAL 
!       ***************************************************

        J=1
        DO 101,K=1,RES
        DATA(J)=P(K,1)
        DATA(J+1)=P(K,2)
        J=J+2
101     CONTINUE
        CALL FOUR1(DATA,NN,1)
        
!       **************************************************
!                DATA IS THE FT OF POTENTIAL. 
!                CALCULATION OF K VECTORS.
!       *************************************************

        DO 25,K=1,RES/2+1 
        W(K)=2.0*PI*REAL(K-1)/REAL(RES)
 25     CONTINUE
        J=1
        DO 26,K=RES/2+2,RES
        W(K)=-1.*W(RES/2+1-J)
        J=J+1
 26     CONTINUE
       
!       ***************************************************
!         FT OF P1 WILL BE CONVERTED IN FT OF ELECTRIC
!         FIELD BY MULTIPLICATION  WITH X (SEE NOTES).
!               PS IS THE FIELD IN THE AIR.
!       ***************************************************

	I=1
        DO 2751,K=1,RES
        R=W(K)
	    PS(I)=R*DATA(I)
	    PS(I+1)=R*DATA(I+1)
	    I=I+2
 2751   CONTINUE

	    I=1
        DO 27,K=1,RES
        R=W(K)
    IF(R.EQ.0.0)THEN
        P1(I)=DATA(I)*(1./D)
        P1(I+1)=DATA(I+1)*(1./D)
    ELSE
		IF(R.GT.0.)THEN
        	P1(I)=DATA(I)*R*(1.+EXP(-2.*R*D))/ &
     	&  (1.-EXP(-2.*R*D))
        	P1(I+1)=DATA(I+1)*R*(1.+EXP(-2.*R*D))/ &
     	&  (1.-EXP(-2.*R*D)) 
		ELSE
	     	P1(I)=DATA(I)*R*(1.+EXP(2.*R*D))/ &
     	&	(EXP(2.*R*D)-1.) 
	     	P1(I+1)=DATA(I+1)*R*(1.+EXP(2.*R*D))/ &
     	&	(EXP(2.*R*D)-1.)
		END IF
    END IF
	    I=I+2
 27     CONTINUE
        
!	*************************************************
!       P1 IS THE FT OF ELECTRIC FIELD.PS IS THE FT OF
!	       THE ELECTRIC FIELD IN THE AIR.
!         INVERSE FT WILL TRANSFORM P1,PS IN ELECTRIC FIELD.
!       **************************************************

	DO 551,I=1,2*RES
	DATA(I)=P1(I)
 551    CONTINUE
        CALL FOUR1(DATA,NN,-1)
	I=1 
	DO 333,K=1,RES
	P1S(K,1)=1./RES*DATA(I)
	P1S(K,2)=1./RES*DATA(I+1)
	I=I+2
 333    CONTINUE

	DO 552,I=1,2*RES
	DATA(I)=PS(I)
 552    CONTINUE
	CALL FOUR1(DATA,NN,-1)
	I=1 
	DO 334,K=1,RES
	PSS(K,1)=1./RES*DATA(I)
	PSS(K,2)=1./RES*DATA(I+1)
	I=I+2
 334    CONTINUE


!       **************************************************
!        CALCULATION OF SUMs(FIELD*SURF.),WHICH IS DESIGNATED
!          BY Cn AND IS THE TOTAL CHARGE ON THE APPROPRIATE 
!        SURFACE DIVIDED BY THE DIELECTRIC CONSTANT OF SILICON.
!          WE HAVE SUMMATION ONLY ON SQUARES  BECAUSE
!              OUTSIDE OF THE SQUARES CHARGE IS ZERO.
!       ***************************************************

        C1=0.
        DO 80,K=A1+1,A2-1
        C1=C1+ESI*P1S(K,1)+EAIR*PSS(K,1)
 80     CONTINUE
        C2=0.
        DO 81,K=A3+1,A4-1
        C2=C2+ESI*P1S(K,1)+EAIR*PSS(K,1)
 81     CONTINUE
        C3=0.
        DO 82,K=A5+1,A6-1
        C3=C3+ESI*P1S(K,1)+EAIR*PSS(K,1)
 82     CONTINUE
        C4=0.    
        DO 83,K=A7+1,A8-1
        C4=C4+ESI*P1S(K,1)+EAIR*PSS(K,1)
 83     CONTINUE
        C5=0.
        DO 84,K=A9+1,A10-1
        C5=C5+ESI*P1S(K,1)+EAIR*PSS(K,1)
 84     CONTINUE

 
!       *****************************************************
!        SETTING ELECTRIC DISPALCEMENT CONTINIOUS IN THE
!        FRONT SURFACE OUTSIDE OF THE JUNCTIONS.
!       *****************************************************

        DO 36,K=1,RES
        IF(((K.LE.A1).OR.((K.GE.A2).AND.(K.LE.A3)) &
     &  .OR.((K.GE.A4).AND.(K.LE.A5)).OR. &
     &  ((K.GE.A6).AND.(K.LE.A7)).OR.((K.GE.A8) &
     &  .AND.(K.LE.A9)).OR.(K.GE.A10))) &
     &  THEN
        P1S(K,1)=EAIR*PSS(K,1)/ESI
        P1S(K,2)=EAIR*PSS(K,2)/ESI
        ELSE
        END IF
 36     CONTINUE
 
!       ***************************************************
!           P1S WILL BE CONVERTED IN FT OF VOLTAGE
!       ***************************************************

	I=1
	DO 335 K=1,RES
	DATA(I)=P1S(K,1)
	DATA(I+1)=P1S(K,2)
	I=I+2
 335    CONTINUE
        CALL FOUR1(DATA,NN,1)
        I=1
        DO 43,K=1,RES
        R=W(K)
        IF(R.EQ.0.)THEN
        DATA(I)=DATA(I)*D
        DATA(I+1)=DATA(I+1)*D
        ELSE
	IF(R.GT.0.)THEN
        DATA(I)=DATA(I)/ &
     &  (R*(1.+EXP(-2.*R*D))/(1.-EXP(-2.*R*D)))
        DATA(I+1)=DATA(I+1)/ &
     &  (R*(1.+EXP(-2.*R*D))/(1.-EXP(-2.*R*D)))
	ELSE
	DATA(I)=DATA(I)/ &
     &	(R*(1.+EXP(2.*R*D))/(EXP(2.*R*D)-1.))
	DATA(I+1)=DATA(I+1)/ &
     &	(R*(1.+EXP(2.*R*D))/(EXP(2.*R*D)-1.))
        END IF
        END IF
	I=I+2
 43     CONTINUE
        
!       ****************************************************
!           DATA IS THE FT OF VOLTAGE.
!        IT WILL BE CONVERTED TO VOLTAGE
!       ****************************************************

        CALL FOUR1(DATA,NN,-1)
	I=1
        DO 234,K=1,RES
        P1S(K,1)=(1./RES)*DATA(I)
        P1S(K,2)=(1./RES)*DATA(I+1)
	I=I+2
 234    CONTINUE
 
!       *******************************************************   
!        NOW P1 IS THE NEW VALUE FOR VOLTAGE.WE SUPPOSE
!        THAT THE SOLUTION IS A SUPERPOSITION OF THE SUPPOSED
!        VALUE ON THE SURFACE AND P1.WE APPLY THE FORMULA:
!            VOLTAGE = ( 1 - f ) P + f P1
!        WHERE f IS A COEFFICIENT BETWEEN 0 AND 1.
!        WE SELECT A VALUE OF f FOR CONVERGENCE OF
!        THE ALGORITHM.
!        NEXT WE CHECK THE DIFFERENCE VOLTAGE(NEW)-VOLTAGE 
!        ALSO WE TAKE CARE OF VOLTAGE TO BE 1 IN THE SQUARE.
!        *******************************************************
 
        DO 52,K=1,RES
        P(K,1)=0.9*P(K,1)+0.1*P1S(K,1)
        P(K,2)=0.9*P(K,2)+0.1*P1S(K,2)
 52     CONTINUE
        DO 193,K=A1+1,A2-1
        P(K,1)=0.
        P(K,2)=0.
 193    CONTINUE
        DO 203,K=A3+1,A4-1
        P(K,1)=0.
        P(K,2)=0.
 203    CONTINUE
        DO 213,K=A5+1,A6-1
        P(K,1)=1.
        P(K,2)=0.
 213    CONTINUE
        DO 223,K=A7+1,A8-1
        P(K,1)=0.
        P(K,2)=0.
 223    CONTINUE
        DO 233,K=A9+1,A10-1
        P(K,1)=0.
        P(K,2)=0. 
 233    CONTINUE

        ICHK=1
        DO 553,K=1,RES
        IF (ABS(P(K,1)-P1S(K,1)) .GT. 1.E-3) ICHK=0
 553    CONTINUE 
       WRITE(22,*)'CYCLE NUMBER ',CYCLE
        C01=(C2+C4)/2.*(-1.)
        C02=(C1+C5)/2.*(-1.)
        C00=C3-2.*C01-2.*C02
        Write(22,*),'Capacitancies, C00, C01, C02'
        WRITE(22,*) C00,C01,C02
	    C1D=LJ*ESI/D
	    C01=C01/C1D
	    C02=C02/C1D
	    C00=C00/C1D
	   Write(22,*),'Capacitancies normalized to C1D=epsilon*(S/D,'
       WRITE(22,*)'C1D=',C1D,'C00=',C00,'C01=',C01,'C02=',C02
       WRITE(22,*)' '
        IF (ICHK .EQ. 0) GO TO 222
        xx=LJ*11.9*8.85e-18/D*1.e4*1.e12
	    print*,C00*xx,C01*xx,C02*xx
       WRITE(22,*)'THE RESOLUTION WAS=',RES,'L=',LJ
       WRITE(22,*)'DEPTH W=',D,'SEPERATION',PITCH
       WRITE(22,*)'RUN No ',Z
       WRITE(22,*)'**********************************************'
 1968   CONTINUE        
        END


	SUBROUTINE FOUR1(DATA,NN,ISIGN)
	DOUBLE PRECISION WR,WI,WPR,WPI,WTEMP,THETA
	DIMENSION DATA(*)
	N=2*NN
	J=1
	DO 11 I=1,N,2
	    IF(J.GT.I)THEN
	       TEMPR=DATA(J)
	       TEMPI=DATA(J+1)
	       DATA(J)=DATA(I)
	       DATA(J+1)=DATA(I+1)
	       DATA(I)=TEMPR
	       DATA(I+1)=TEMPI
	    ENDIF
	    M=N/2
  1         IF ((M.GE.2).AND.(J.GT.M)) THEN
		J=J-M
		M=M/2
            GO TO 1
	    ENDIF
	    J=J+M
 11 	    CONTINUE
        MMAX=2
  2	IF (N.GT.MMAX) THEN
	    ISTEP=2*MMAX
	    THETA=6.28318530717959D0/(ISIGN*MMAX)
	    WPR=-2.D0*DSIN(0.5D0*THETA)**2
	    WPI=DSIN(THETA)
            WR=1.D0
	    WI=0.D0
	    DO 13 M=1,MMAX,2
	    DO 12 I=M,N,ISTEP
		J=I+MMAX
		TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
                TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J) 
		DATA(J)=DATA(I)-TEMPR
                DATA(J+1)=DATA(I+1)-TEMPI
		DATA(I)=DATA(I)+TEMPR
		DATA(I+1)=DATA(I+1)+TEMPI
 12		 CONTINUE
	    WTEMP=WR
	    WR=WR*WPR-WI*WPI+WR
	    WI=WI*WPR+WTEMP*WPI+WI
 13	     CONTINUE
	 MMAX=ISTEP
	GO TO 2
	ENDIF
	RETURN
	END
