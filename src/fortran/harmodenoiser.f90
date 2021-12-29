module harmodenoiser
! ------------------------------------------------------------------------------
!
! ------------------------------------------------------------------------------
use calculus
contains
! ------------------------------------------------------------------------------
subroutine hd_nbnt_(nb,nt_,nt__,nt,fo,dt)
  ! ----------------------------------------------------------------------------
  ! nt : total number of time samples
  ! fo : central frequncy to look for         ( nt__ = ceil(1/fo/dt) )
  ! dt : dt in time
  ! nb : desired # of time blocks
  !
  ! example for nt=4000, fo=50, dt=2.5e-4, nb=4
  !        [nb,nt_,nt__] = hd_nbnt_(4000,50,2.5e-4,4)
  ! nb=4, nt_=180, nt__=80.
  ! ----------------------------------------------------------------------------
  ! dt_ : desired time duration of time-blocks (dt_=(nt*dt)/(desired # of blocks))
  ! nt_ : number of samples in the over-lapping time
  ! ----------------------------------------------------------------------------
  ! nb  : # of time blocks
  ! nt_ : number of samples in the over-lapping time
  ! nt__: overlapping number of time samples   (nt__ = ceil(1/fo/dt))
  ! ----------------------------------------------------------------------------
  !
  ! the forward model is,
  !
  ! uh = cos_blocs · α + sin_blocs · β
  !
  ! the cos_blocs matrix (nt × nb·nh) looks a little like this:
  !
  !                       nh
  !                       ↓
  !                    _____________________________
  !        |          |        |                    |
  !        |  nt_ →   |    *   |_________     0     |
  !        |          |________| ← nt__  |          |
  ! nt →   |          |        |    *    |          | · α
  !        |          |        |_________|          |
  !        |          |   0                  etc    |
  !        |          |_____________________________|
  !
  !                                  ↑
  !                                nb·nh
  !
  ! and one cos_bloc looks like:
  !
  !             nh
  !             ↓
  !          ________
  !         |        |
  ! nt_ →   |    *   | = cos( 2*pi*fo*t*h )
  !         |________|
  !
  ! each block is of size nt_ x nh.
  ! they all overlap on nt__ samples.
  ! this big matrix is of size nt x nb*nh.
  ! there are nb = (nt-nt__)/(nt_-nt__) blocks.
  ! α is of size nb*nh x 1.
  ! ----------------------------------------------------------------------------
  ! --
  ! nt_  :: size(1) of blocks
  ! --
  ! nt__ :: size(1) of overlap
  ! take dt * nt__ = 1/fo
  !
  !     --> nt__ = 1/fo/dt
  ! --
  ! nb   :: # of blocks
  ! each block is of size nt_ x nh
  ! nt = (# of blocks)*(size(1) of blocks) - collisions
  ! nt = nb*nt_ - nt__*(nb - 1)
  !
  ! -->    nb = (nt-nt__)/(nt_-nt__)
  !
  ! nb should always be an integer!
  !
  ! assuming nt__ is fixed (because nt__ depends on the frequency to be found),
  ! we can change nt_ (always an integer) until we get an integer for nb.
  ! we have,
  !
  ! nt_ = ((nt-nt__) / nb) + nt__
  !
  ! so nb must divide nt-nt__.
  !
  ! we can factor nt-nt__ and then look for an nb which satisfies our estimate
  ! length for dt_ = nt_*dt,
  !
  ! dt_ ≈ (((nt-nt__)/(some combination of factor(nt-nt__) )) + nt__) * dt
  !
  ! nb = the best combination of factor(nt-nt__) that is close to dt_.
  !
  ! whatever that combination may be, it may change the initial value for nt_!
  !
  ! nt_ = ((nt-nt__) / nb) + nt__.
  !
  ! ----------------------------------------------------------------------------
  integer, intent(in) :: nt
  double precision, intent(in) :: fo, dt
  integer, intent(in out) :: nb, nt_, nt__

  integer, allocatable :: fact_combi(:), ierr_(:)
  integer :: ifact, ifact_, nfact
  double precision, allocatable :: err_(:)
  ! ----------------------------------------------------------------------------
  ! overlapping number of time samples
  nt__ = ceiling(1/fo/dt)
  ! get number of factors of ntnt
  ntnt  = nt-nt__
  nfact = 0
  do ifact = 1,floor(dsqrt(dble(ntnt)))
    if (mod(ntnt,ifact)==0) then
      if (ntnt/ifact == ifact) then
        nfact = nfact+1
      else
        nfact = nfact+2
      endif
    endif
  enddo
  ! ----------------------------------------------------------------------------
  ! get list of these factors
  allocate(fact_combi(nfact))
  ifact_ = 1
  do ifact = 1,floor(dsqrt(dble(ntnt)))
    if (mod(ntnt,ifact)==0) then
      if (ntnt/ifact == ifact) then
        fact_combi(ifact_) = ifact
        ifact_=ifact_+1
      else
        fact_combi(ifact_) = ifact
        fact_combi(ifact_+1) = ntnt/ifact
        ifact_=ifact_+2
      endif
    endif
  enddo
  ! ----------------------------------------------------------------------------
  ! get nb by matching to dt_
  ! desired time duration of time-blocks (dt_=(nt*dt)/(desired # of blocks))
  dt_ = (nt*dt)/nb
  ! error bucket
  allocate(err_(nfact))
  do ifact=1,nfact
    ! estimate dt_
    ! dt_esti_ = (((nt-nt__)/( some combination of factor(nt-nt__) )) + nt__) * dt
    dt_esti_ = dble( ((nt-nt__)/( fact_combi(ifact) )) + nt__ ) * dt
    err_(ifact) = abs(dt_esti_ - dt_)
  enddo
  allocate(ierr_(nfact))
  do ifact=1,nfact
    ierr_(ifact) = ifact
  enddo
  call quicksort(err_,ierr_,nfact)
  ifact = ierr_(1) ! min index
  nb = fact_combi(ifact)
  ! ----------------------------------------------------------------------------
  ! number of samples in the over-lapping time
  nt_ = ((nt-nt__) / nb) + nt__
end subroutine hd_nbnt_
! ------------------------------------------------------------------------------
subroutine hd_fwd(uh,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__
  double precision, intent(in) :: t(nt)
  double precision, intent(in) :: alphas(nb*nh), betas(nb*nh), fos(nb), h(nh)
  double precision, intent(in out) :: uh(nt)

  integer :: it_, ih, indexes_h(nh), indexes_t(nt_)
  double precision :: argu, cos_bloc(nt_,nh), sin_bloc(nt_,nh), uh_(nt)
  double precision, parameter :: pi=3.14159265358979
  ! ----------------------------------------------------------------------------
  do ib=1,nb
    ! --- one block ---
    ! each block is of size nt_ × nh
    ! cos_bloc = cos( 2*pi*fo*t*h )

    ! time-interval times
    ! indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
    do it_=1,nt_
      indexes_t(it_) = it_ + (ib-1)*(nt_-nt__)
    enddo
    ! harmonic interval
    ! indexes_h = (1 + (ib-1)*nh):(nh + (ib-1)*nh);
    do ih=1,nh
      indexes_h(ih) = ih + (ib-1)*nh
    enddo

    ! each bloc is a matrix of size nt_ × nh
    do it_=1,nt_
      do ih=1,nh
        argu = 2*pi*fos(ib)*t(indexes_t(it_)) * h(ih)
        cos_bloc(it_,ih) = dcos(argu)
        sin_bloc(it_,ih) = dsin(argu)
      enddo
    enddo

    do it_=1,nt
      uh_(it_) = 0
    enddo

    ! uh_(indexes_t) = cos_bloc * alphas(indexes_h) + sin_bloc * betas(indexes_h);
    do it_=1,nt_
      argu = 0
      do ih=1,nh
        argu = argu + cos_bloc(it_,ih) * alphas(indexes_h(ih)) + sin_bloc(it_,ih) * betas(indexes_h(ih))
      enddo
      uh_(indexes_t(it_)) = argu
    enddo

    ! record the harmonic signal WITH overlapping blocks
    ! uh = uh + uh_;
    do it_=1,nt
      uh(it_) = uh(it_) + uh_(it_)
    enddo
  enddo
end subroutine hd_fwd
! ------------------------------------------------------------------------------
subroutine hd_grad_f()
  
  ! ----------------------------------------------------------------------------


end subroutine hd_grad_f
! ------------------------------------------------------------------------------
end module harmodenoiser
