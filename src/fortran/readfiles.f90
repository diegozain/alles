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
subroutine read_dbl(dbl,filename)
  character(*), intent(in) :: filename
  double precision, intent(in out) :: dbl

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)dbl
  close(1)
end subroutine read_dbl
! ------------------------------------------------------------------------------
subroutine read_vec_sgl(mati,nlen,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen
  real, intent(in out) :: mati(nlen)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1)mati
  close(1)
end subroutine read_vec_sgl
! ------------------------------------------------------------------------------
subroutine read_mat_sgl(mati,nrows,ncols,filename)
  character(*), intent(in) :: filename
  integer, intent(in) :: nrows,ncols
  real, intent(in out) :: mati(nrows,ncols)

  open(unit=1,file=filename,access='stream',form='unformatted')
  read(1) mati
  close(1)
end subroutine read_mat_sgl
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
subroutine save_vec_int(mati,nlen,filename,id)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen, id
  integer, intent(in out) :: mati(nlen)
  integer :: imati

  open(unit=id,file=filename,access='stream',form='unformatted')
  write(id) (mati(imati),imati=1,nlen)
  close(id)
end subroutine save_vec_int
! ------------------------------------------------------------------------------
subroutine save_vec_dbl(mati,nlen,filename,id)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen, id
  double precision, intent(in out) :: mati(nlen)
  integer :: imati

  open(unit=id,file=filename,access='stream',form='unformatted')
  write(id) (mati(imati),imati=1,nlen)
  close(id)
end subroutine save_vec_dbl
! ------------------------------------------------------------------------------
subroutine save_vec_cmplx(mati,nlen,filename,id)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlen, id
  complex*16, intent(in out) :: mati(nlen)
  integer :: imati

  open(unit=id,file=filename,access='stream',form='unformatted')
  write(id) (mati(imati),imati=1,nlen)
  close(id)
end subroutine save_vec_cmplx
! ------------------------------------------------------------------------------
subroutine read_filey(filey,nlines,filename,id)
  character(*), intent(in) :: filename
  integer, intent(in) :: nlines, id
  character*256, intent(in out) :: filey(nlines)

  character*256 :: chartmp
  integer :: iline = 0

  open(unit=id,file=filename)
  do iline = 1, nlines
    read(id,'(A)') filey(iline)
  end do

  close(id)
end subroutine read_filey
! ------------------------------------------------------------------------------
subroutine get_nfiley(nlines,filename,id)
  character(*), intent(in) :: filename
  integer, intent(in) :: id
  integer, intent(in out) :: nlines

  character*256 :: chartmp
  integer :: ierr = 0

  nlines = 0
  open(unit=id,file=filename)
  do while (ierr == 0)
    nlines = nlines + 1
    read(id,*,iostat=ierr) chartmp
  end do
  nlines = nlines - 1

  close(id)
end subroutine get_nfiley
! ------------------------------------------------------------------------------
end module readfiles
