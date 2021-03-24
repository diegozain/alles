module fourier
! ------------------------------------------------------------------------------
! diego domenzain @ CSM 2021
!
! this module is about fourier.
!
! * 'dfork' is a discrete fast fourier transform by Thomas Forbriger.
! * 'costap' is a cosine taper by Thomas Forbriger.
!
! Thomas Forbriger's code: https://git.scc.kit.edu/Seitosh/Seitosh
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! the following fast fourier transform is just copied from
! the program refseis from j. ungerer
!
! this code was originally published by gerhard müller
! in his lecture notes on digital signal processing:
! gerhard mueller, 1999. digitale signalverarbeitung. skriptum zur
! gleichnamigen vorlesung. institut für meteorologie und geophysik,
! universität frankfurt.
!
! the original algorithm appears to be due to claerbout, j.f., 
! "fundamentals of geophysical data processing with applications 
! to petroleum prospecting", mcgraw-hill, 1976.
!
!
! 23456789012345678901234567890123456789012345678901234567890123456789012
! *** subroutine fuer fouriertransformation ************************** !
! die von gerherd mueller verwendetet schnelle fouriertransformation   !
! fork wurde umgeschrieben fuer double complex                         !
! es muessen implementiert sein: double complex,dcmplx,cdexp           !
!                                                                      !
! zum verfahren der schnellen fouriertransformation(fft) und zur ar-   !
! beitsweise von fork siehe g.mueller: digitale signalverarbeitung i,  !
! vorlesungsmanuskript.                                                !
!                                                                      !
! variablen:                                                           !
!    lx       seismogrammlaenge bzw. anzahl der stuetzstellen,abtast-  !
!             werte des seismogramms/spektrums.muss eine zeier-potenz  !
!             sein.                                                    !
!    cx(lx)   feld auf dessen realteil die funktionswerte der zeit-    !
!             funktion stehen und nach transformation ihre fourierko-  !
!             effizienten.                                             !
!    signi    signi=-1.d0 bedeutet berechnung der fourierkoeffizienten !
!             signi=+1.d0 bedeutet ruecktransformation                 !
! ******************************************************************** !     

      subroutine dfork(lx,cx,signi)
      integer     i,istep,j,l,lx,m
      real*8      sc,pi,signi
      complex*16  cx(lx),carg,cw,ctemp

      pi=3.14159265358979d0
      j=1
      sc=1.d0/dble(lx)
      sc=dsqrt(sc)
      do 5  i=1,lx
      if(i.gt.j) goto 2
      ctemp=cx(j)*sc
      cx(j)=cx(i)*sc
      cx(i)=ctemp
2     m=lx/2
3     if(j.le.m) goto 5
      j=j-m
      m=m/2
      if(m.ge.1) goto 3
5     j=j+m
      l=1
6     istep=2*l
      do 8  m=1,l
      carg=dcmplx(0.,1.)*(pi*signi*dble(m-1))/dble(l)
      cw=cdexp(carg)
      do 8 i=m,lx,istep
      ctemp=cw*cx(i+l)
      cx(i+l)=cx(i)-ctemp
8     cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.lx) goto 6
      return
      end subroutine dfork
! ------------------------------------------------------------------------------
double precision function costap(min,wil,wir,max,val)
double precision min,wil,wir,max,val
double precision pi
parameter(pi=3.1415926535898D0)
if (val.lt.min) then
 costap=0.D0
elseif(val.le.wil) then
 if (wil.eq.min) then
    costap=0.D0
 else
    costap=0.5D0-0.5D0*dcos((val-min)*pi/(wil-min))
 endif
elseif(val.gt.max) then
 costap=0.D0
elseif(val.ge.wir) then
 if (max.eq.wir) then
    costap=0.D0
 else
    costap=0.5D0+0.5D0+dcos((val-wir)*pi/(max-wir))
 endif
else
 costap=1.D0
endif
return
end function costap
! ------------------------------------------------------------------------------
end module fourier