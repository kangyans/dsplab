%% PSF simulation for different spiral trajectory
%
%  Note:
% --------------
%   This is only the simplest case :-( where no signal decay was considered
%


clear; close all;clc;

addpath(genpath('D:\UVA\ReconSimu\MRIRecon_Toolbox\vds'));
addpath(genpath('D:\UVA\ReconSimu\MRIRecon_Toolbox\Spiral_deblurring_Toolbox\matlab-programs\SpiralDeblur'));

Nint    = 4        ;
fov     = 220      ;
res     = 1        ;
imsize  = [220,220];

RR =  [0,1,1/8];

for i = 1:length(RR)

    [k{i},w{i}]   = retraced_spiral_generator(Nint,fov,res,RR(i));
    psf(:,:,i)    = nufft_recon(ones(size(k{i})),k{i},w{i},imsize);
    mpt(i)        = max(abs(psf(imsize(1)/2+1,:,i)));
    profile(:,i)  = abs(psf(imsize(1)/2+1,:,i))./mpt(i);
end

sp1 = k{1,1};
sp2 = k{1,2};
sp3 = k{1,3};

nx  = size(sp1,1)     ;
ptn = ceil(RR(3)*nx) ;

figure(1),
subplot(131),plotf(real(sp1(:,1)),imag(sp1(:,1)),'r--');
subplot(132),plotf(real(sp2(1:nx,1)),imag(sp2(1:nx,1)),'b-');hold on;
plotf(real(sp2(nx:end,1)),imag(sp2(nx:end,1)),'r--');hold off;
subplot(133),plotf(real(sp3(1:ptn,1)),imag(sp3(1:ptn,1)),'b-');hold on;
plotf(real(sp3(ptn:end,1)),imag(sp3(ptn:end,1)),'r--');hold on;


figure(2),
plot(profile1,'b-' ,'LineWidth',3);
hold on;
plot(profile2,'g-' ,'LineWidth',3);
plot(profile2,'r-' ,'LineWidth',3);
hold off;
legend('Spiral out', 'Spiral in-out','Asymetric Spiral in-out');
xlabel('# of points');
ylabel('Point Spread Function (a.u.)');
set(gca,'FontSize',20,'FontWeight','bold');
axis square;


rmpath('D:\UVA\ReconSimu\MRIRecon_Toolbox\Spiral_deblurring_Toolbox\matlab-programs\SpiralDeblur');
rmpath('D:\UVA\ReconSimu\MRIRecon_Toolbox\vds');


