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
 integer :: it, ib, ih, ibh
 integer, parameter :: nhyper=11
 double precision :: hyperparam(nhyper)

 ! 💾
 character(len=128) :: path_save, name_

 ! ⌚
 real :: start_time, end_time, rate_time
 ! -----------------------------------------------------------------------------
 path_save = 'bin/'
 ! -----------------------------------------------------------------------------
 dt = 2.5e-4 ! sec
 fo = 5      ! Hz
 nt = 8000   ! # of samples 8000

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
 !
 !                                 🎵
 !
 ! -----------------------------------------------------------------------------
 ! ⌚
 start_time = omp_get_wtime()
 ! -----------------------------------------------------------------------------
 call harmodenoi(uo,dt,fo,h,nt,nh)
 ! -----------------------------------------------------------------------------
 ! ⌚
 end_time = omp_get_wtime()
 print *, 'elapsed time: ', (end_time - start_time), 'seconds'
 ! -----------------------------------------------------------------------------
 print*,''
 print*,'this is slower wtf???'
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
 ! -----------------------------------------------------------------------------
 !                            😎😎 bye bye ✋✋
 ! -----------------------------------------------------------------------------
 print *,''
 print *,' --- bye bye --- '
 print *,''
end program harmodenoi_synt
