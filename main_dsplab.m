%  
%  (d)iffustion-weighted (SP)iral imaging simulation (LAB). Ispired by
%  the work of Feizollah, S et al. see the reference below:
%  (www.github.com/TardifLab/dMRI_sequence_simulations)
%
%  Kang Yan (2023), UVa.
%


clear 
clc
close all


% -------------------------------------------------------------------------
% Add toolbox
% -------------------------------------------------------------------------
addpath(genpath('D:\UVA\ReconSimu\srs\sp\nufft_recon'));
addpath(genpath('D:\UVA\ReconSimu\srs\sp\spiralDeblur'));
addpath(genpath('D:\UVA\ReconSimu\srs\sim\vds'));
addpath(genpath('D:\UVA\ReconSimu\srs\sim\epg'));


% -------------------------------------------------------------------------
% set sequence parameters
% -------------------------------------------------------------------------

seq.B0           = 3    ;  % main field magnitude [unit: Tesla]
seq.TR           = 3000 ;  % time repitition [unit: ms]
seq.TE           = 50   ;  % echo time       [unit: ms]
seq.fov          = 0.256;  % tield of View   [unit: m]
seq.G_max        = 40e-3;  % maximum diffusion-encoding gradient [unit: T/m]
seq.dRF_exc      = 2.56 ;  % RF excitation duration [unit: ms]
seq.dRF_refoc    = 2.56 ;  % RF refocsing duration  [unit: ms]
seq.num_dummy    = 4    ;  % number of dummy sequences


% -------------------------------------------------------------------------
% Set simulation parameters
% -------------------------------------------------------------------------
seq.decaymodel   = 'epg'   ;
seq.res          = 1e-3    ;   % image resolution [unit: m]
seq.tissue       = 'WM'    ;   % the tissue in the brain

seq.R            = 1       ;   % acceleration ratio
seq.PF           = 1       ;   % partial fourier 
seq.BW           = 1384    ;   % bandwidth (Hz/px)
seq.b_value      = 1000*1e6;   % b-value
seq.BWpp         = 1e6     ;
seq.extraR       = 1/6     ;
seq.Nint         = 1       ;
seq.output       = 'image' ;
seq.addnoise     = 0       ;
N  = floor(seq.fov./seq.res);
seq.imsize = [N,N];

% -------------------------------------------------------------------------
% Set off resonance parameters
% -------------------------------------------------------------------------
seq.offres_mode   = 1           ;
seq.offres_fmap   = off_res_fmap_sim(seq.imsize, 100);
seq.offres_method = 'lstseg'    ;
seq.offres_nseg   = 10          ;
seq.B0correction  = 1           ;
seq.B0correction_method = 'icpr' ;

% -------------------------------------------------------------------------
% Set imaging object parameters
% -------------------------------------------------------------------------
if strcmp(seq.output,'psf')
    seq.image = zeros(N,N);
    seq.image(end/2+1,end/2+1) = 1;
end

if strcmp(seq.output,'image')    
    [GM, WM, CSF, mask] = phantom_generator(N) ;
    im  = csphantom(N); 
    im  = phantom(N);
    seq.image = im;
end

% -------------------------------------------------------------------------
% Simulation
% -------------------------------------------------------------------------
seq.traj_type = 'sp';
sp   = seq;
sp   = diff_spiral_signal_model(sp);

seq.traj_type = 'spoi';
spoi   = seq;
spoi   = diff_spiral_signal_model(spoi);

seq.traj_type = 'spaio';
spaio = seq;
spaio.TE = 140;
spaio = diff_spiral_signal_model(spaio) ;


% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
if strcmp(seq.output, 'psf') == 1
    figure(10),
    
    subplot(121),
    plot(abs(sp.psf(:,end/2+1)),'b'); hold on;
    plot(abs(spoi.psf(:,end/2+1)),'r');
    plot(abs(spaio.psf(:,end/2+1)),'g');hold off;
    legend('Spiral out', 'Spiral out-in','Spiral in-out-in');
    title('Point Spread Function');
    axis square
    
    subplot(122),
    plot(1:length(sp.kdata),abs(sp.kdata(:,1)),'b');hold on;
    plot(1:length(spoi.kdata),abs(spoi.kdata(:,1)),'r');
    plot(1:length(spaio.kdata),abs(spaio.kdata(:,1)),'g');
    hold off;
    legend('Spiral out', 'Spiral out-in','Spiral in-out-in');
    title('Signal decay along spiral readout');
    axis square
end

if strcmp(seq.output,'image') == 1 && (seq.offres_mode == 1)
    figure(20),imshow(abs([sp.psf spoi.psf spaio.psf; sp.psf2 spoi.psf2 spaio.psf2]),[]);
elseif strcmp(seq.output,'image') == 1 && (seq.offres_mode == 0)
    figure(20),imshow(abs([sp.psf spoi.psf spaio.psf]),[]);
end
