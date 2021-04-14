program get_data
 ! -----------------------------------------------------------------------------
 ! diego domenzain @ CSM 2021
 ! 
 ! generate simple data, and then have python plot it.
 ! -----------------------------------------------------------------------------
 implicit none
 ! -----------------------------------------------------------------------------
 double precision, parameter :: pi=3.1415926535898
 complex*16, parameter :: IMAGINE=(0.D0,1.D0)
 
 ! 1D vector
 integer, parameter :: nt = 25600
 double precision :: dt = 0.001D0
 double precision t(nt), one_d(nt)
 integer :: it

 ! 2D matrix
 integer, parameter :: nrows=3, ncols=2
 integer :: irows,icols,ientry
 integer two_d(nrows,ncols)
 ! -----------------------------------------------------------------------------
 ! 1D vector
 ! true source init 
 do it=1,nt
   t(it)=dble(it-1)*dt
   one_d(it)=dsin(t(it))
 end do
 ! -----------------------------------------------------------------------------
 ! ! save as characters
 ! open(1, file = 't.dat', status='unknown')
 ! open(2, file = 'one_d.dat', status='unknown')
 ! do it = 1,nt
 !  write(1,"(E15.7)") t(it)
 !  write(2,"(E15.7)") one_d(it)
 ! end do  
 ! close(1)
 ! close(2)
 ! -----------------------------------------------------------------------------
 ! save as binary
 open(1,file = 't.dat', status='unknown',form='unformatted',access='stream')
 open(2,file = 'one_d.dat', status='unknown',form='unformatted',access='stream')
 write(1) (t(it),it=1,nt)
 write(2) (one_d(it),it=1,nt)
 close(1)
 close(2)
 ! -----------------------------------------------------------------------------
 ! 2D matrix 
 ientry=1
 do icols=1,ncols
   do irows=1,nrows
     two_d(irows,icols)=ientry
     ientry=ientry+1
   end do
 end do
 ! -----------------------------------------------------------------------------
 ! ! save as characters
 ! open(1, file = 'two_d.dat', status='unknown')
 ! do icols=1,ncols
 !   do irows=1,nrows
 !     write(1,"(I2)") two_d(irows,icols)
 !   end do
 ! end do
 ! close(1)
 ! -----------------------------------------------------------------------------
 ! save as binary
 open(1,file = 'two_d.dat', status='unknown',form='unformatted',access='stream')
 write(1) ((two_d(irows,icols),irows=1,nrows),icols=1,ncols)
 close(1)
 ! -----------------------------------------------------------------------------
end program get_data