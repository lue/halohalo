c**********************************************************************
c          SOME USEFUL AUXILIARY FUNCTIONS & SUBROUTINES
c
c   - subroutine for linear interpolation
c   - ARCTANH
c   - TANH
c   - XEXP
c   - lblnk
c
c  Frank van den Bosch                                    Dec. 1999
c**********************************************************************

      SUBROUTINE linintpol(xa,ya,N,x,y)

      IMPLICIT NONE

      INTEGER N,j
      REAL    xa(N),ya(N),x,y,x1,x2,y1,y2

c---

      CALL locate(xa,N,x,j)
      x1 = xa(j)
      x2 = xa(j+1)
      y1 = ya(j)
      y2 = ya(j+1)

      y = y1 + ((x-x1)/(x2-x1)) * (y2-y1)

      RETURN
      END      

c**********************************************************************

      REAL FUNCTION arctanh(x)
c--------------------------------------------------------------------
c
c Auxialiary function to compute arctanh(x)
c
c--------------------------------------------------------------------          
        
      REAL    x

c---

      arctanh = 0.5 * ALOG((1.0+x)/(1.0-x))

      END

c**********************************************************************

      REAL FUNCTION TANH(x)
c--------------------------------------------------------------------
c
c Auxialiary function to compute TANH(x)
c
c--------------------------------------------------------------------

      REAL     x

      REAL     XEXP
      EXTERNAL XEXP

c---

      TANH = (XEXP(x) - XEXP(-x)) / (XEXP(x) + XEXP(-x)) 

      END

c**********************************************************************

      REAL FUNCTION XEXP(x)
c--------------------------------------------------------------------
c
c Auxialiary function to compute EXP(x)
c
c--------------------------------------------------------------------

      REAL    x

c---

      IF (x.LT.-40.0) THEN
        XEXP = 0.0
      ELSE
        XEXP = EXP(x)
      END IF

      END

c**********************************************************************

      INTEGER FUNCTION lblnk(char)
c--------------------------------------------------------------------
c
c Function gives NUMBER of characters in a character variable `char'
c
c--------------------------------------------------------------------

      IMPLICIT NONE

      character char*(*)

c---

      lblnk=index(char,' ')-1

      RETURN
      END

c**********************************************************************

      SUBROUTINE Terminate(message)
c--------------------------------------------------------------------
c
c  Output error message and terminate program
c
c--------------------------------------------------------------------

      IMPLICIT NONE

      character  message*(*)

c---

      WRITE(*,'(A)')message

      STOP

      RETURN
      END

c**********************************************************************

      SUBROUTINE conf_levels(yy,n,s02,s16,s50,s84,s98)
c---------------------------------------------------------------------------
c
c  Given an array yy(1:n), sort and determine the 
c  5, 16, 50, 84 and 95 percent confidence levels.
c
c---------------------------------------------------------------------------

      INTEGER n,n02,n16,n50,n84,n98
      REAL    yy(n),xn,s02,s16,s50,s84,s98

c---

      CALL sort(n,yy)

      xn = FLOAT(n)       
      n02 = NINT(0.0228 * xn)
      n16 = NINT(0.1587 * xn)
      n50 = NINT(0.5 * xn)
      n84 = NINT(0.8413 * xn)
      n98 = NINT(0.9772 * xn)

      s02 = yy(n02)
      s16 = yy(n16)
      s50 = yy(n50)
      s84 = yy(n84)
      s98 = yy(n98)

      RETURN
      END
      
c**********************************************************************

      SUBROUTINE conf_levels2(yy,n,s33,s67,n33,n67)
c---------------------------------------------------------------------------
c
c  Given an array yy(1:n), sort and determine the 
c  33 and 67 percent confidence levels.
c
c---------------------------------------------------------------------------

      INTEGER n,n33,n67
      REAL    yy(n),xn,s33,s67

c---

      CALL sort(n,yy)

      xn = FLOAT(n) 
            
      n33 = NINT(0.3333 * xn)
      n67 = NINT(0.6667 * xn)

      s33 = yy(n33)
      s67 = yy(n67)

      RETURN
      END

c**********************************************************************






