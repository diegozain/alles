module readfiles
! ------------------------------------------------------------------------------
! diego domenzain 2021
!
! read binary files.
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine size_mat(mati_size,dime,filename)
  character(*), intent(in) :: filename
  integer, intent(in out) :: mati_size(3), dime

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)mati_size
  close(1)

  dime = 0
  do idime=1,3
    dime = dime+1
    if (mati_size(idime) .eq. 0) then
      dime = dime-1
      exit
    endif
  enddo
end subroutine size_mat
! ------------------------------------------------------------------------------
subroutine read_mat_dbl(mati,nrows,ncols,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nrows,ncols
  double precision, intent(in out) :: mati(nrows,ncols)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1) mati
  close(1)
end subroutine read_mat_dbl
! ------------------------------------------------------------------------------
subroutine read_vec_dbl(mati,nlen,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen
  double precision, intent(in out) :: mati(nlen)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)mati
  close(1)
end subroutine read_vec_dbl
! ------------------------------------------------------------------------------
subroutine read_mat_int(mati,nrows,ncols,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nrows,ncols
  integer, intent(in out) :: mati(nrows,ncols)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)mati
  close(1)
end subroutine read_mat_int
! ------------------------------------------------------------------------------
subroutine read_vec_int(mati,nlen,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen
  integer, intent(in out) :: mati(nlen)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)mati
  close(1)
end subroutine read_vec_int
! ------------------------------------------------------------------------------
subroutine save_vec_dbl(mati,nlen,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen
  double precision, intent(in out) :: mati(nlen)
  integer :: imati

  open(unit=1,file=filename,access='stream',form='unformatted')
  write(1) (mati(imati),imati=1,nlen)
  close(1)
end subroutine save_vec_dbl
! ------------------------------------------------------------------------------
end module readfiles
