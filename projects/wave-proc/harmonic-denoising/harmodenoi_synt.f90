program harmodenoi_synt
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
 use omp_lib, only : omp_get_wtime
 implicit none
 ! -----------------------------------------------------------------------------
 ! 🎵
 integer :: nt,nb,nh,nt_,nt__,nw
 double precision :: dt, fo
 double precision, allocatable :: uo(:), h(:), alphas(:), betas(:), fos(:), t(:)
 double precision, allocatable :: us(:)

 ! 📟
 integer :: it, ib, ih, ibh, nhyper
 double precision, allocatable :: hyperparam(:)

 ! 💾
 character(len=128) :: path_save, name_

 ! ⌚
 real :: start_time, end_time, rate_time
 ! -----------------------------------------------------------------------------
 path_save = 'bin\'
 ! -----------------------------------------------------------------------------
 dt = 2.5e-4 ! sec
 fo = 5      ! Hz
 nt = 8000   ! # of samples

 nh = 4
 nb = 1
 ! -----------------------------------------------------------------------------
 ! input  :: target nb, fo, & dt
 ! output :: new nb, nt_, & nt__
 call hd_nbnt_(nb,nt_,nt__,nt,fo,dt)
 ! -----------------------------------------------------------------------------
 print*,''
 print*,'----------------------------------------------------------------------'
 print*,'      total time     = ',nt*dt,' seconds'
 print*,'      interval time  = ',nt_*dt,' seconds'
 print*,'      overlap time   = ',nt__*dt,' seconds'
 print*,''
 print*,'      # of blocks nb = ',nb
 print*,'----------------------------------------------------------------------'
 print*,''
 print*,''
 ! -----------------------------------------------------------------------------
 allocate(uo(nt))
 allocate(h(nh))
 allocate(t(nt))

 allocate(alphas(nb*nh))
 allocate(betas(nb*nh))
 allocate(fos(nb))

 allocate(us(nt))
 ! -----------------------------------------------------------------------------
 ! build time
 do it=1,nt
   t(it) = (it-1)*dt
 enddo
 ! -----------------------------------------------------------------------------
 !                         👉 synthetic fwd 👉
 ! -----------------------------------------------------------------------------
 alphas= [2,1,1,1,2,1,1,1]
 betas = [2,1,1,1,2,1,1,1]
 h     = [1,3,5,11]
 do ib=1,nb
   fos(ib) = fo
 enddo
 ! do some harmonic
 call hd_fwd(uo,t,alphas,betas,fos,h,nt,nb,nh,nt_,nt__)
 ! put in some hidden signal
 do it=1,nt
   us(it) = dexp(-t(it)**2)
   uo(it) = uo(it) + us(it)
 enddo
 ! -----------------------------------------------------------------------------
 ! 💾 uo & t
 name_ = trim(path_save) // 'uo_obs.bin'
 call save_vec_dbl(uo,nt,name_)
 name_ = trim(path_save) // 'uo_sig.bin'
 call save_vec_dbl(us,nt,name_)
 ! -----------------------------------------------------------------------------
 !                          📟 hyperparam 📟
 ! -----------------------------------------------------------------------------
 nhyper = 11
 allocate(hyperparam(nhyper))
 ! -----------------------------------------------------------------------------
 ! k_fos_ & k_fos__
 hyperparam(1) = 1e-8
 hyperparam(2) = 2e-5
 ! k_alphas_ & k_alphas__
 hyperparam(3) = 1e-8
 hyperparam(4) = 1e-3
 ! k_betas_ & k_betas__
 hyperparam(5) = 1e-8
 hyperparam(6) = 5e-3
 ! nparabo_fos & nparabo_a & nparabo_b
 hyperparam(7) = 105
 hyperparam(8) = 105
 hyperparam(9) = 105
 ! niter_fos & niter_ab
 hyperparam(10) = 5
 hyperparam(11) = 5
 ! -----------------------------------------------------------------------------
 !                              fo stuff
 ! -----------------------------------------------------------------------------
 ! window to convolve with
 nw = ceiling((1/(fo*h(1)))/dt)

 ! fo varies very little between blocks.
 ! in fact, with EM data i am yet to see a case where nb > 1 gives better
 ! results than nb=1.
 do ib=1,nb
   fos(ib) = fo
 enddo
 ! -----------------------------------------------------------------------------
 !                                 🎵
 ! -----------------------------------------------------------------------------
 ! ⌚
 start_time = omp_get_wtime()
 ! -----------------------------------------------------------------------------
 call harmodenoi_(uo,t,h,alphas,betas,fos,nt,nb,nh,nt_,nt__,nw,&
   hyperparam, nhyper)
 ! -----------------------------------------------------------------------------
 ! ⌚
 end_time = omp_get_wtime()
 print *, 'elapsed time: ', (end_time - start_time), 'seconds'
 ! -----------------------------------------------------------------------------
 !👣
 print*,''
 print*,'--> recovered frequencies per block:'
 print*,''
 do ib=1,nb
   print*,fos(ib),'(Hz)'
 enddo

 print*,''
 print*,'--> recovered alphas per block:'
 print*,''
 ibh=1
 do ib=1,nb
   do ih=1,nh
     print*,alphas(ibh),'(alpha)'
     ibh = ibh+1
   enddo
   print*,''
 enddo

 print*,''
 print*,'--> recovered betas per block:'
 print*,''
 ibh=1
 do ib=1,nb
   do ih=1,nh
     print*,betas(ibh),'(beta)'
     ibh = ibh+1
   enddo
   print*,''
 enddo
 print*,''
 ! -----------------------------------------------------------------------------
 ! 💾 uo & t
 name_ = trim(path_save) // 'uo_rec.bin'
 call save_vec_dbl(uo,nt,name_)
 name_ = trim(path_save) // 't.bin'
 call save_vec_dbl(t,nt,name_)
 ! -----------------------------------------------------------------------------
 deallocate(uo)
 deallocate(h)
 deallocate(t)
 deallocate(us)

 deallocate(alphas)
 deallocate(betas)
 deallocate(fos)

 deallocate(hyperparam)
 ! -----------------------------------------------------------------------------
 !                            😎😎 bye bye ✋✋
 ! -----------------------------------------------------------------------------
 print *,''
 print *,' --- bye bye --- '
 print *,''
end program harmodenoi_synt
