module harmodenoiser
! ------------------------------------------------------------------------------
! diego domenzain
!
! 🔥 super 🍧 cool 🎵 harmonic denoising.
!
! based on my own matlab code,
! which is in turn based on a seismic processing paper from Geophysics.
! -- put name of that paper here --
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
  !        [nb,nt_,nt__] = hd_nbnt_(nt,fo,dt,nb)
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
  if (nb==1) then
    nt__ = 0
  else
    nt__ = ceiling(1/(fo*dt))
  endif
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

  integer :: ib, it_, ih, index_h, index_t
  double precision :: argu, dotter, uh_(nt)
  double precision, parameter :: pi=3.14159265358979
  ! ----------------------------------------------------------------------------
  do it_=1,nt
    uh(it_)=0
  enddo
  do ib=1,nb
    ! --- one block ---
    ! each block is of size nt_ × nh
    ! cos_bloc = cos( 2*pi*fo*t*h )

    ! initialize uh_
    do it_=1,nt
      uh_(it_) = 0
    enddo

    ! each bloc is a matrix of size nt_ × nh
    do it_=1,nt_
      dotter = 0
      index_t = it_ + (ib-1)*(nt_-nt__)
      do ih=1,nh
        argu   = 2*pi*fos(ib)*t(index_t) * h(ih)
        ! uh_ = cos_bloc·α + sin_bloc·β
        index_h = ih + (ib-1)*nh
        dotter = dotter + dcos(argu)*alphas(index_h) + dsin(argu)*betas(index_h)
      enddo
      uh_(index_t) = dotter
    enddo

    ! record the harmonic signal WITH overlapping blocks
    ! uh = uh + uh_;
    do it_=1,nt
      uh(it_) = uh(it_) + uh_(it_)
    enddo
  enddo
end subroutine hd_fwd
! ------------------------------------------------------------------------------
subroutine hd_grad_f(g_fos,error_,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__
  double precision, intent(in) :: error_(nt), t(nt)
  double precision, intent(in) :: alphas(nb*nh), betas(nb*nh), fos(nb), h(nh)
  double precision, intent(in out) :: g_fos(nb)

  integer :: ib, it_, ih, index_h, index_t
  double precision :: argu, argu_, dotter, del_fo(nt_)
  double precision, parameter :: pi=3.14159265358979
  ! ----------------------------------------------------------------------------
  do ib=1,nb
    ! each bloc is a matrix of size nt_ × nh
    do it_=1,nt_
      dotter = 0
      index_t = it_ + (ib-1)*(nt_-nt__)
      do ih=1,nh
        argu =  2*pi*fos(ib)*t(index_t) * h(ih)
        argu_=  2*pi*t(index_t) * h(ih)
        index_h = ih + (ib-1)*nh
        dotter = dotter - argu_*dsin(argu)*alphas(index_h) + argu_*dcos(argu)*betas(index_h)
      enddo
      del_fo(it_) = dotter
    enddo

    ! dot product
    ! g_fos(ib) = error_(indexes_t) * (cos_bloc__ + sin_bloc__);
    dotter = 0
    do it_=1,nt_
      index_t= it_ + (ib-1)*(nt_-nt__)
      dotter = dotter + error_(index_t) * del_fo(it_)
    enddo
    g_fos(ib) = dotter
  enddo
end subroutine hd_grad_f
! ------------------------------------------------------------------------------
subroutine hd_grad_a(g_alphas,error_,t,fos,h,nt,nb,nh,nt_,nt__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__
  double precision, intent(in) :: error_(nt), t(nt)
  double precision, intent(in) :: fos(nb), h(nh)
  double precision, intent(in out) :: g_alphas(nb*nh)

  integer :: ib, ibh, it_, ih, index_t, nbh
  double precision :: argu, cos_bloc_, dotter
  double precision, parameter :: pi=3.14159265358979
  ! ----------------------------------------------------------------------------
  nbh = nb*nh
  do ibh=1,nbh
   ! translate one for-loop into two for-loops
   ih = mod(ibh,nh)
   if (ih==0) then
     ih=nh
   endif
   ib = ((ibh - ih) / nh) + 1

   ! time-interval times
   ! indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
   ! &
   ! cosine block
   ! bloc is a matrix of size nt_ × 1
   ! &
   ! dot product
   dotter = 0
   do it_=1,nt_
     index_t = it_ + (ib-1)*(nt_-nt__)

     argu = 2*pi*fos(ib)*t(index_t) * h(ih)
     cos_bloc_ = dcos(argu)

     dotter = dotter + error_(index_t) * cos_bloc_
   enddo
   g_alphas(ibh) = dotter
 enddo
end subroutine hd_grad_a
! ------------------------------------------------------------------------------
subroutine hd_grad_b(g_betas,error_,t,fos,h,nt,nb,nh,nt_,nt__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__
  double precision, intent(in) :: error_(nt), t(nt)
  double precision, intent(in) :: fos(nb), h(nh)
  double precision, intent(in out) :: g_betas(nb*nh)

  integer :: ib, ibh, it_, ih, index_t, nbh
  double precision :: argu, sin_bloc_, dotter
  double precision, parameter :: pi=3.14159265358979
  ! ----------------------------------------------------------------------------
  nbh = nb*nh
  do ibh=1,nbh
   ! translate one for-loop into two for-loops
   ih = mod(ibh,nh)
   if (ih==0) then
     ih=nh
   endif
   ib = ((ibh - ih) / nh) + 1

   ! time-interval times
   ! indexes_t = (1 + (ib-1)*(nt_-nt__)):(nt_+ (ib-1)*(nt_-nt__));
   ! &
   ! sine block
   ! bloc is a matrix of size nt_ × 1
   ! &
   ! dot product
   dotter = 0
   do it_=1,nt_
     index_t = it_ + (ib-1)*(nt_-nt__)

     argu = 2*pi*fos(ib)*t(index_t) * h(ih)
     sin_bloc_ = dsin(argu)

     dotter = dotter + error_(index_t) * sin_bloc_
   enddo
   g_betas(ibh) = dotter
 enddo
end subroutine hd_grad_b
! ------------------------------------------------------------------------------
subroutine hd_obj(objfnc_,error_,data_, datao_, OBJ, nd)
  integer, intent(in) :: nd
  double precision, intent(in) :: data_(nd), datao_(nd)
  integer, intent(in) :: OBJ
  double precision, intent(in out) :: objfnc_, error_(nd)
  ! ----------------------------------------------------------------------------
  ! OBJ =
  !       0 : sum of squared residuals        Θ = ∑ e
  !       1 : log( sum of squared residuals ) Θ = ∑ log(e)
  !       2 : entropy                         Θ = - ∑ e^2 ⋅ log( e^2 )
  ! ----------------------------------------------------------------------------
  objfnc_ = 0
  do id_=1,nd
    error_(id_) = data_(id_) - datao_(id_)
    objfnc_ = objfnc_ + error_(id_)**2
  enddo

  if (OBJ == 1) then
    objfnc_ = dlog(objfnc_)
    do id_=1,nd
      error_(id_) = (1/objfnc_) * error_(id_)
    enddo
  elseif (OBJ == 2) then
    objfnc_ = 0
    do id_=1,nd
      objfnc_ = objfnc_ - (error_(id_)**2) * dlog(error_(id_)**2)
    enddo
    do id_=1,nd
      error_(id_) = - error_(id_) - error_(id_)*dlog(error_(id_)**2)
    enddo
  endif
end subroutine hd_obj
! ------------------------------------------------------------------------------
subroutine hd_step_f(step_fos,uo,g_fos,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
  nparabo,OBJ,k_fos_,k_fos__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__, nparabo, OBJ
  double precision, intent(in) :: uo(nt), t(nt), g_fos(nb)
  double precision, intent(in) :: alphas(nb*nh), betas(nb*nh), fos(nb), h(nh)
  double precision, intent(in) :: k_fos_,k_fos__
  double precision, intent(in out) :: step_fos

  double precision :: Ob_(nparabo), k_fos(nparabo), fos_(nb), Ob
  double precision :: uh_(nt), error_(nt)
  integer :: iparabo, ib, iOb, iOb_(nparabo)
  ! ----------------------------------------------------------------------------
  ! --- build perturbations
  k_fos = logspace(dlog10(k_fos_),dlog10(k_fos__),nparabo)
  ! compute many objective function values
  do iparabo=1,nparabo
    ! perturb
    ! fo_ = fo ⊙ exp(-k·g_fo ⊙ fo)
    do ib=1,nb
      fos_(ib) = fos(ib)*dexp(- k_fos(iparabo)*g_fos(ib)*fos(ib))
    enddo
    ! fwd
    call hd_fwd(uh_,t,alphas,betas,fos_,h,nt,nb,nh,nt_,nt__)
    ! obj
    call hd_obj(Ob, error_, uh_, uo, OBJ, nt)
    Ob_(iparabo) = Ob
  enddo
  ! -- brute line-search
  do iOb=1,nparabo
    iOb_(iOb) = iOb
  enddo
  call quicksort(Ob_,iOb_,nparabo)
  iOb = iOb_(1) ! min index
  step_fos = k_fos(iOb)
end subroutine hd_step_f
! ------------------------------------------------------------------------------
subroutine hd_step_a(step_alphas,uo,g_alphas,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
  nparabo,OBJ,k_alphas_,k_alphas__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__, nparabo, OBJ
  double precision, intent(in) :: uo(nt), t(nt), g_alphas(nb*nh)
  double precision, intent(in) :: alphas(nb*nh), betas(nb*nh), fos(nb), h(nh)
  double precision, intent(in) :: k_alphas_,k_alphas__
  double precision, intent(in out) :: step_alphas

  double precision :: Ob_(nparabo), k_alphas(nparabo), alphas_(nb*nh), Ob
  double precision :: uh_(nt), error_(nt)
  integer :: iparabo, ib, iOb, iOb_(nparabo)
  ! ----------------------------------------------------------------------------
  ! --- build perturbations
  k_alphas = linspace(k_alphas_,k_alphas__,nparabo)
  ! compute many objective function values
  do iparabo=1,nparabo
    ! perturb
    ! α_ = α - k·g_α
    do ib=1,nb*nh
      alphas_(ib) = alphas(ib) - k_alphas(iparabo)*g_alphas(ib)
    enddo
    ! fwd
    call hd_fwd(uh_,t,alphas_,betas,fos,h,nt,nb,nh,nt_,nt__)
    ! obj
    call hd_obj(Ob, error_, uh_, uo, OBJ, nt)
    Ob_(iparabo) = Ob
  enddo
  ! -- brute line-search
  do iOb=1,nparabo
    iOb_(iOb) = iOb
  enddo
  ! !🐛
  ! do ib=1,nparabo
  !   print*,Ob_(ib),iOb_(ib),k_alphas(ib)
  ! enddo
  call quicksort(Ob_,iOb_,nparabo)
  iOb = iOb_(2) ! min index
  ! !🐛
  ! do ib=1,nparabo
  !   print*,Ob_(ib),iOb_(ib)
  ! enddo
  step_alphas = k_alphas(iOb)
end subroutine hd_step_a
! ------------------------------------------------------------------------------
subroutine hd_step_b(step_betas,uo,g_betas,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
  nparabo,OBJ,k_betas_,k_betas__)
  integer, intent(in) :: nt, nb, nh, nt_, nt__, nparabo, OBJ
  double precision, intent(in) :: uo(nt), t(nt), g_betas(nb*nh)
  double precision, intent(in) :: alphas(nb*nh), betas(nb*nh), fos(nb), h(nh)
  double precision, intent(in) :: k_betas_,k_betas__
  double precision, intent(in out) :: step_betas

  double precision :: Ob_(nparabo), k_betas(nparabo), betas_(nb*nh), Ob
  double precision :: uh_(nt), error_(nt)
  integer :: iparabo, ib, iOb, iOb_(nparabo)
  ! ----------------------------------------------------------------------------
  ! --- build perturbations
  k_betas = linspace(k_betas_,k_betas__,nparabo)
  ! compute many objective function values
  do iparabo=1,nparabo
    ! perturb
    ! β_ = β - k·g_β
    do ib=1,nb*nh
      betas_(ib) = betas(ib) - k_betas(iparabo)*g_betas(ib)
    enddo
    ! fwd
    call hd_fwd(uh_,t,alphas,betas_,fos,h,nt,nb,nh,nt_,nt__)
    ! obj
    call hd_obj(Ob, error_, uh_, uo, OBJ, nt)
    Ob_(iparabo) = Ob
  enddo
  ! -- brute line-search
  do iOb=1,nparabo
    iOb_(iOb) = iOb
  enddo
  call quicksort(Ob_,iOb_,nparabo)
  iOb = iOb_(1) ! min index
  step_betas = k_betas(iOb)
end subroutine hd_step_b
! ------------------------------------------------------------------------------
subroutine hd_hyperparam(k_fos_,k_fos__,k_alphas_,k_alphas__,k_betas_,k_betas__,&
  nparabo_fos,nparabo_a,nparabo_b,niter_fos,niter_ab,hyperparam,nhyper)
  integer, intent(in) :: nhyper
  double precision, intent(in) :: hyperparam(nhyper)
  integer, intent(in out) :: nparabo_fos,nparabo_a,nparabo_b,niter_fos,niter_ab
  double precision, intent(in out) :: k_fos_,k_fos__
  double precision, intent(in out) :: k_alphas_,k_alphas__,k_betas_,k_betas__
  ! ----------------------------------------------------------------------------
  k_fos_     = hyperparam(1)
  k_fos__    = hyperparam(2)
  k_alphas_  = hyperparam(3)
  k_alphas__ = hyperparam(4)
  k_betas_   = hyperparam(5)
  k_betas__  = hyperparam(6)
  nparabo_fos= int(hyperparam(7))
  nparabo_a  = int(hyperparam(8))
  nparabo_b  = int(hyperparam(9))
  niter_fos  = int(hyperparam(10))
  niter_ab   = int(hyperparam(11))
end subroutine hd_hyperparam
! ------------------------------------------------------------------------------
subroutine harmodenoi_(uo,t,h,alphas,betas,fos,nt,nb,nh,nt_,nt__,nw,&
  hyperparam, nhyper)
  integer, intent(in) :: nt, nb, nh, nt_, nt__, nw
  double precision, intent(in) :: t(nt), h(nh), hyperparam(nhyper)
  double precision, intent(in out) :: uo(nt),alphas(nh*nb),betas(nh*nb),fos(nb)

  integer :: ibh, ib, ih, iter
  double precision :: dt, x_, objfnc_, error_(nt)

  double precision :: k_fos_,k_fos__,k_alphas_,k_alphas__,k_betas_,k_betas__
  integer :: nparabo_fos,nparabo_a,nparabo_b,niter_fos,niter_ab

  double precision, allocatable :: fos_niter(:,:), a_niter(:,:), b_niter(:,:)
  double precision, allocatable :: ob_fos(:), ob_ab(:)
  integer, allocatable :: iob_fos(:), iob_ab(:)

  double precision :: uh(nt), g_alphas(nh*nb), g_betas(nh*nb), g_fos(nb)
  double precision :: step_alphas, step_betas, step_fos
  ! ----------------------------------------------------------------------------
  !                    💠 initial guess for α & β & fos 💠
  ! ----------------------------------------------------------------------------
  ! the idea for this one is that:
  ! • α & β are the ones contributing most of the noise in the signal (std)
  ! • α & β for small frequencies are usually larger when dealing with EM data (/ibh)
  call std(x_,uo,nt)
  ibh=1
  do ib=1,nb
    do ih=1,nh
      alphas(ibh)= x_ / ih
      betas(ibh) = x_ / ih
      ibh = ibh+1
    enddo
  enddo
  ! ----------------------------------------------------------------------------
  !                            📟 hyperparam 📟
  ! ----------------------------------------------------------------------------
  call hd_hyperparam(k_fos_,k_fos__,k_alphas_,k_alphas__,k_betas_,k_betas__,&
    nparabo_fos,nparabo_a,nparabo_b,niter_fos,niter_ab,hyperparam,nhyper)
  ! ok now lets allocate some fail-safe stuff
  allocate(fos_niter(nb,niter_fos))
  allocate(a_niter(nb*nh,niter_ab))
  allocate(b_niter(nb*nh,niter_ab))
  allocate(ob_fos(niter_fos))
  allocate(ob_ab(niter_ab))
  allocate(iob_fos(niter_fos))
  allocate(iob_ab(niter_ab))
  ! ----------------------------------------------------------------------------
  !
  !                      💃🎸 lets rock and roll 🎸💃
  !
  ! ----------------------------------------------------------------------------
  !                                🎹 fo 🎹
  ! ----------------------------------------------------------------------------
  do iter=1,niter_fos
    ! fwd & obj
    call hd_fwd(uh,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
    call hd_obj(objfnc_,error_,uh, uo, 1, nt)
    ! fail-safe
    ob_fos(iter) = objfnc_
    do ib=1,nb
      fos_niter(ib,iter) = fos(ib)
    enddo
    ! gradient
    call hd_grad_f(g_fos,error_,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
    ! step-size
    call hd_step_f(step_fos,uo,g_fos,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
      nparabo_fos,1,k_fos_,k_fos__)
    ! update (always positive)
    do ib=1,nb
      fos(ib) = fos(ib) * dexp(-step_fos*g_fos(ib)*fos(ib))
    enddo
  enddo
  ! init this guy 👲
  do ib=1,niter_fos
    iob_fos(ib) = ib
  enddo
  call quicksort(ob_fos,iob_fos,niter_fos)
  iter = iob_fos(1) ! min index
  do ib=1,nb
    fos(ib) = fos_niter(ib,iter)
  enddo
  ! ----------------------------------------------------------------------------
  !                                 α & β
  ! ----------------------------------------------------------------------------
  do iter=1,niter_ab
    ! fwd & obj
    call hd_fwd(uh,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
    call hd_obj(objfnc_,error_,uh, uo, 0, nt)
    ! fail-safe
    ob_ab(iter) = objfnc_
    do ibh=1,nb*nh
      a_niter(ibh,iter) = alphas(ibh)
      b_niter(ibh,iter) = betas(ibh)
    enddo
    ! α gradient
    call hd_grad_a(g_alphas,error_,t,fos,h,nt,nb,nh,nt_,nt__)
    ! α step-size
    call hd_step_a(step_alphas,uo,g_alphas,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
      nparabo_a,0,k_alphas_,k_alphas__)
    ! β gradient
    call hd_grad_b(g_betas,error_,t,fos,h,nt,nb,nh,nt_,nt__)
    ! β step-size
    call hd_step_b(step_betas,uo,g_betas,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__,&
      nparabo_b,0,k_betas_,k_betas__)
    ! update
    do ibh=1,nb*nh
      alphas(ibh)= alphas(ibh)- step_alphas*g_alphas(ibh)
      betas(ibh) = betas(ibh) - step_betas*g_betas(ibh)
    enddo
  enddo
  ! init this guy 👲
  do ibh=1,niter_ab
    iob_ab(ibh) = ibh
  enddo
  call quicksort(ob_ab,iob_ab,niter_ab)
  iter = iob_ab(1) ! min index
  do ibh=1,nb*nh
    alphas(ibh)= a_niter(ibh,iter)
    betas(ibh) = b_niter(ibh,iter)
  enddo
  ! ----------------------------------------------------------------------------
  !                          👉 last pass 👉
  ! ----------------------------------------------------------------------------
  call hd_fwd(uh,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
  ! overwrite uo with solution
  do it=1,nt
    uo(it) = uo(it) - uh(it)
  enddo
  call window_mean(uo,nt,nw)
end subroutine harmodenoi_
! ------------------------------------------------------------------------------
end module harmodenoiser
