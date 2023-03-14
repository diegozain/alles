program harmodenoibin
 ! -----------------------------------------------------------------------------
 ! diego domenzain 2021
 !
 !
 !                     🎵 harmonic denoising demo 🎵
 !
 !
 ! -----------------------------------------------------------------------------
 use harmodenoiser
 use readfiles
 use omp_lib
 implicit none
 ! -----------------------------------------------------------------------------
 ! 🎵
 integer :: nt, nb, nh, nt_, nt__, nw, nabmn, nttail, nttrim
 double precision :: dt, fo
 double precision, allocatable :: h(:), alphas(:), betas(:), fos(:), t(:), uo(:)
 real, allocatable :: dataips(:,:)
 double precision, allocatable :: alphas_(:), betas_(:), fos_(:)
 real, allocatable :: dataips_(:)

 ! 📟
 integer :: it, ih, ib, ibh, iabmn, nn
 integer, parameter :: nhyper=11
 double precision :: hyperparam(nhyper)

 ! 💾 👓
 character*256, allocatable :: filey(:)
 character(len=128) :: path_save, path_read, path_read_, name_
 integer :: dime, mati_size(3), nlines

 ! ⬜⬛ openMP
 ! in 💩 & 🚀 uncomment this one, comment next1
 integer, parameter :: nthreads=10
 ! integer :: nthreads
 integer :: stacksize

 ! ⌚
 real :: start_time, end_time
 integer :: date_time(8)
 character*10 date_time_(3)
 ! -----------------------------------------------------------------------------
 call date_and_time(date_time_(1), date_time_(2), date_time_(3), date_time)
 print *,''
 print *,'                            🎵🎵🎵🎵'
 print *,''
 print *, '--------------------------------------------------------------------'
 print *,''
 print '("  ",i2,".",i2,".",i4)',date_time(3),date_time(2),date_time(1)
 print '("  ",i2,":",i2,":",i2)',date_time(5),date_time(6),date_time(7)
 ! -----------------------------------------------------------------------------
 !
 !                                   👁️ 👃 👁️
 !
 ! -----------------------------------------------------------------------------
 call get_nfiley(nlines,'paths.txt',2)
 allocate(filey(nlines))
 call read_filey(filey,nlines,'paths.txt',2)
 path_save = filey(1)
 path_read = filey(2)
 deallocate(filey)

 print *,''
 print*,'----------------------------------------------------------------------'
 print *,''
 print *, '📝  reading from:   ', path_read
 print *, '💾  saving in:      ', path_save
 ! -----------------------------------------------------------------------------
 path_read_ = trim(path_read) // 'dataips_size.bin'
 call size_mat(mati_size,dime,path_read_)
 nt = mati_size(1)
 nabmn = mati_size(2)
 allocate(dataips(nt,nabmn))
 path_read_ = trim(path_read) // 'dataips.bin'
 call read_mat_sgl(dataips,nt,nabmn,path_read_)
 ! -----------------------------------------------------------------------------
 !                          📟 hyperparam 📟
 ! -----------------------------------------------------------------------------
 ! k_fos_ & k_fos__
 hyperparam(1) = 1e-9
 hyperparam(2) = 1e-4
 ! k_alphas_ & k_alphas__
 hyperparam(3) = 1e-8
 hyperparam(4) = 1e-2
 ! k_betas_ & k_betas__
 hyperparam(5) = 1e-8
 hyperparam(6) = 1e-2
 ! nparabo_fos & nparabo_a & nparabo_b
 hyperparam(7) = 20
 hyperparam(8) = 50
 hyperparam(9) = 50
 ! niter_fos & niter_ab
 hyperparam(10) = 6
 hyperparam(11) = 6
 ! -----------------------------------------------------------------------------
 dt = 2.5e-4 ! sec
 fo = 9      ! Hz

 nh = 6
 nb = 1
 ! -----------------------------------------------------------------------------
 !
 !                                   👷👷👷👷
 !
 ! -----------------------------------------------------------------------------
 allocate(t(nt))
 ! build ⏰
 do it=1,nt
   t(it) = (it-1)*dt
 enddo
 ! 🎼
 allocate(h(nh))
 do ih=1,nh
   h(ih) = ih
 enddo
 ! 🏠 window to convolve with
 nw = ceiling((1/(fo*h(1)))/dt)
 ! ✂️ ⏰ ✂️
 nttail = ceiling(0.1*nh*(1/fo)/dt)
 nttrim = nt-nttail
 ! -----------------------------------------------------------------------------
 ! input  :: target nb, fo, & dt
 ! output :: new nb, nt_, & nt__
 call hd_nbnt_(nb,nt_,nt__,nt,fo,dt)

 allocate(dataips_(nttrim*nabmn))
 allocate(alphas_(nb*nh*nabmn))
 allocate(betas_(nb*nh*nabmn))
 allocate(fos_(nb*nabmn))
 ! -----------------------------------------------------------------------------
 print*,'----------------------------------------------------------------------'
 print*,'      number of time-series is ',nabmn
 print*,'      each with                ',nt,' samples'
 print*,'----------------------------------------------------------------------'
 print*,'      total time     = ',nt*dt,' seconds'
 print*,'      interval time  = ',nt_*dt,' seconds'
 print*,'      overlap time   = ',nt__*dt,' seconds'
 print*,''
 print*,'         # of blocks = ',nb
 print*,'      # of harmonics = ',nh
 print*,'            harmonic = ',fo,' Hz'
 print*,'----------------------------------------------------------------------'
 print*,''
 ! -----------------------------------------------------------------------------
 !                    ⬜⬛◻️◼️ openMP ◻️◼️⬜⬛
 ! -----------------------------------------------------------------------------
 ! in 💩 & 🚀 uncomment this one, comment next1
 call omp_set_num_threads(nthreads)
!  nthreads = omp_get_max_threads()
 call kmp_set_stacksize_s(500000) ! 8000000 = 8 Gb
 stacksize = kmp_get_stacksize_s()
 print *, '    # of threads is',nthreads
 print *, '    stack size = ',dble(stacksize)/10**6,'Gb'
 print*,''
 ! -----------------------------------------------------------------------------
 !
 !                                 🎵🎵🎵🎵
 !
 ! -----------------------------------------------------------------------------
 ! ⌚
 start_time = omp_get_wtime()
 ! -----------------------------------------------------------------------------
 !$omp parallel private(iabmn,it,ib,ibh) &
 !$omp firstprivate(uo,alphas,betas,fos) &
 !$omp shared(dataips_,alphas_,betas_,fos_) &
 !$omp shared(dataips,t,h,hyperparam) &
 !$omp shared(nabmn,nt,nttrim,nb,nh,nt_,nt__,nw)
 allocate(uo(nt))
 allocate(alphas(nb*nh))
 allocate(betas(nb*nh))
 allocate(fos(nb))
 !$omp do schedule(static)
 do iabmn=1,nabmn
   do it=1,nt
     uo(it) = dble(dataips(it,iabmn))
   enddo
   do ib=1,nb
     fos(ib) = fo
   enddo
   ! 🎵⛔
   call harmodenoi_(uo,t,h,alphas,betas,fos,nt,nb,nh,nt_,nt__,nw,hyperparam,nhyper)
   ! 🔽🔀
   do it=1,nttrim
     dataips_(it + (iabmn-1)*nttrim) = real(uo(it))
   enddo
   do ibh=1,nb*nh
     alphas_(ibh + (iabmn-1)*nb*nh) = alphas(ibh)
     betas_(ibh + (iabmn-1)*nb*nh) = betas(ibh)
   enddo
   do ib=1,nb
     fos_(ib + (iabmn-1)*nb) = fos(ib)
   enddo
 enddo
 !$omp end do
 deallocate(uo)
 deallocate(alphas)
 deallocate(betas)
 deallocate(fos)
 !$omp end parallel
 ! -----------------------------------------------------------------------------
 ! 💾
 name_ = trim(path_save) // 'dataips__size.bin'
 mati_size(1) = nttrim
 mati_size(2) = nabmn
 call save_vec_int(mati_size,3,name_,1)
 name_ = trim(path_save) // 'dataips_.bin'
 call save_vec_sgl(dataips_,nttrim*nabmn,name_,1)
 name_ = trim(path_save) // 'bafos_size.bin'
 mati_size(1) = nb
 mati_size(2) = nh
 mati_size(3) = nabmn
 call save_vec_int(mati_size,3,name_,1)
 name_ = trim(path_save) // 'alphas_.bin'
 call save_vec_dbl(alphas_,nb*nh*nabmn,name_,1)
 name_ = trim(path_save) // 'betas_.bin'
 call save_vec_dbl(betas_,nb*nh*nabmn,name_,1)
 name_ = trim(path_save) // 'fos_.bin'
 call save_vec_dbl(fos_,nb*nabmn,name_,1)
 ! -----------------------------------------------------------------------------
 ! 🚿
 deallocate(t)

 deallocate(h)
 deallocate(dataips)

 deallocate(dataips_)
 deallocate(alphas_)
 deallocate(betas_)
 deallocate(fos_)
 ! -----------------------------------------------------------------------------
 ! ⌚
 end_time = omp_get_wtime()
 print *, 'elapsed time: ', (end_time - start_time), 'seconds'
 ! -----------------------------------------------------------------------------
 !                            😎😎 bye bye ✋✋
 ! -----------------------------------------------------------------------------
 print *,''
 print *,' --- bye bye --- '
 print *,''
end program harmodenoibin
