#include "/Applications/MATLAB_R2014b.app/extern/include/fintrf.h"
C======================================================================
#if 0
C     
C     hhd_solve.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C diego domenzain
C @ Colorado School of Mines
C   May 2021
C
C     hhd_solve.f
C
C     solve for helmholtz-hodge decomposition.
C
C            in Matlab this will be:
C                 [phi,psi] = hhd_solve(x,y,udiv,urot);
C
C I based it off of this example:
C     /Applications/MATLAB_R2014b.app/extern/examples/refbook/xtimesy.F
C      
C     This is a MEX-file for MATLAB.
C     
C======================================================================
C     Gateway routine
C      nlhs: # of input arguments
C      plhs: pointer to output arguments
C      nrhs: # of output arguments
C      prhs: pointer to input arguments
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetN

C     Pointers to input/output mxArrays:
      mwPointer x_ptr, y_ptr, udiv_ptr, urot_ptr, phi_ptr, psi_ptr

C     Array information:
      mwPointer nx, ny
      mwSize nxy

C     Arguments for computational routine:
      real*8 x(1000), y(1000), udiv(250000), urot(250000)
      real*8 phi(250000), psi(250000)

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 4) then
         call mexErrMsgIdAndTxt('MATLAB:hhd_solve:nInput',
     +                           'four inputs required.')
      elseif(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt('MATLAB:hhd_solve:nOutput',
     +                           'Too many output arguments.')
      endif

C     Get the size of the input array.
C
C              n = mxGetN( * )
C           .-----------------.
C           |                 |
C           |                 | m = mxGetM( * )
C           |                 |
C           .-----------------.
C
      nx = mxGetN(prhs(1))
      ny = mxGetN(prhs(2))
      nxy= nx*ny

C     Create Fortran array from the input argument.
      x_ptr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(x_ptr,x,nx)
      y_ptr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(y_ptr,y,ny)
      udiv_ptr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(udiv_ptr,udiv,nxy)
      urot_ptr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(urot_ptr,urot,nxy)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(ny,nx,0)
      phi_ptr = mxGetPr(plhs(1))
      plhs(2) = mxCreateDoubleMatrix(ny,nx,0)
      psi_ptr = mxGetPr(plhs(2))

C     Call the computational subroutine.
      call hhd_solve(nx,ny,x,y,udiv,urot,phi,psi)

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(phi,phi_ptr,nxy)
      call mxCopyReal8ToPtr(psi,psi_ptr,nxy)

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine hhd_solve(nx,ny,x,y,udiv,urot,phi,psi)
C            in Matlab this will be:
C                 [phi,psi] = hhd_solve(x,y,urot,udiv);
      implicit none
      
      real*8, parameter :: pi=3.14159265358979
      
      integer*8, intent(in)  :: nx, ny
      real*8, intent(in)  :: x(nx), y(ny)
      real*8, intent(in)  :: udiv(ny,nx), urot(ny,nx)
      real*8, intent(out) :: phi(ny,nx), psi(ny,nx)
      
      integer*8 ix, iy, ixo, iyo
      real*8 g, phi_, psi_, dx, dy
      
      dx = x(2)-x(1)
      dy = y(2)-y(1)
      
      do ixo=1,nx
            do iyo=1,ny
                  phi_=0.0
                  psi_=0.0
                  do ix=1,nx
                        do iy=1,ny
                              
C                             take care of singularity
                              if (ixo==ix .and. iyo==iy) then
                                    g = 0.0
                              else
      g=(-1.0/(2.0/pi))*dlog(dsqrt((x(ix)-x(ixo))**2+(y(iy)-y(iyo))**2))
                              endif
                              
C                             integral
                              phi_=phi_ + g*udiv(iy,ix)*dx*dy
                              psi_=psi_ + g*urot(iy,ix)*dx*dy
                        enddo
                  enddo
                  phi(iyo,ixo) = phi_
                  psi(iyo,ixo) = psi_
            enddo
      enddo
      return
      end
