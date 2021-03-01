dt=2.95e-03;

load('vz_f')                      
figure;fancy_imagesc(vz_(:,:,700));title('source f')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

vz_=differentiate_cube(vz_,dt);
figure;fancy_imagesc(vz_(:,:,700));title('diff(source f)')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

load('vz_fdt')                                            
figure;fancy_imagesc(vz_(:,:,700));title('source dt(f)')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

vz_=integrate_cube(vz_,dt);                               
figure;fancy_imagesc(vz_(:,:,700));title('integral(source dt(f))')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

