function params=seq_sim_main_v2(params)

% >>>>>>>>>> define relaxation times <<<<<<<<<<
%
% T1, T2, T2* in (ms) used from:
%
%   [1].Cox and Gowland , 2010
%   [2].Peters    et al., 2007
%   [3].Rooney    et al., 2007
%   [4].Wansapura et al., 1999

params.T1_7T_wm  = 1200;
params.T1_7T_gm  = 2000;
params.T1_3T_wm  = 800;
params.T1_3T_gm  = 1300;
params.T2_3T_wm  = 79.6;
params.T2_3T_gm  = 72;
params.T2_7T_wm  = 47;
params.T2_7T_gm  = 47;
params.T2s_7T_wm = 26.8;
params.T2s_7T_gm = 33.2;
params.T2s_3T_gm = 66;
params.T2s_3T_wm = 53.2;

% >>>>>>>>>> set relaxation times based on tissue and field strength <<<<<<<<<<

if(params.B0==3 && params.tissue == "WM")
    params.T1  = params.T1_3T_wm ;
    params.T2  = params.T2_3T_wm ;
    params.T2s = params.T2s_3T_wm;
elseif(params.B0==3 && params.tissue == "GM")
    params.T1  = params.T1_3T_gm;
    params.T2  = params.T2_3T_gm;
    params.T2s = params.T2s_3T_gm;
elseif(params.B0==7 && params.tissue == "WM")
    params.T1  = params.T1_7T_wm;
    params.T2  = params.T2_7T_wm;
    params.T2s = params.T2s_7T_wm;
elseif(params.B0==7 && params.tissue == "GM")
    params.T1  = params.T1_7T_gm;
    params.T2  = params.T2_7T_gm;
    params.T2s = params.T2s_7T_gm;
else
    error('relaxation parameters not defined....')
end

gamma = 267.522e6;    % gyromagnetic ratio in (rad/s/T)

mat_size = ceil(params.fov/params.res);  % nominal matrix size
BW       = params.BWpp*mat_size       ;  % total bandwidth

% >>>>>>>>>> generate spiral trajectory <<<<<<<<<<

if(params.traj_type == "sp")
        
    params.kloc = seq_sim_spiral_generator_v2(params.Nint,params.res,params.fov,params.traj_type,params.extraR);
    
    % >>>>>>>>>> calculate diffusion-encoding duration <<<<<<<<<<
    
    C = [2/3*gamma^2*params.G_max^2,gamma^2*params.G_max^2*(params.RF_refoc_dur*1e-3),0,-params.b_value];
    root         = roots(C);
    params.tdiff = root(find(real(root)>0))*1000;
    pre_readout  = ceil(params.tdiff + params.RF_refoc_dur/2) + ceil(params.RF_exc_dur/2);
    
    params.TE   = 2*pre_readout;
    sample_time = 1e-3; % 1 us 
    params.time_acquisition = params.TE + (0:sample_time:length(params.kloc)*sample_time - sample_time)';
    
    params.readout_dur = length(params.kloc)*sample_time;
    params.label       = strcat('Spiral');
    params.kloc_PF     = [];


elseif(params.traj_type == 'spaio')

    [params.kloc,~,~,params.extraPtn] = seq_sim_spiral_generator_v2(params.Nint, params.res, params.fov,params.traj_type,params.extraR);

    C = [2/3*gamma^2*params.G_max^2,gamma^2*params.G_max^2*(params.RF_refoc_dur*1e-3),0,-params.b_value];
    root         = roots(C);
    params.tdiff = root(find(real(root)>0))*1000;
    sample_time  = 1e-3;
    
    params.time_acquisition = params.TE + (-params.extraPtn*sample_time:sample_time:(length(params.kloc)-params.extraPtn)*sample_time - sample_time)';
    
    params.readout_dur = length(params.kloc)*sample_time;
    params.label       = strcat('Spiral Asymmetric in-out');
    params.kloc_PF     = [];
    
% >>>>>>>>>> generate EPI trajectory <<<<<<<<<<
    
elseif(params.traj_type=="ep")
    
    [params.kloc,params.time_acquisition,dur_to_TE,~,params.kloc_PF,~]=...
        seq_sim_EPI_generator(params.fov,ceil(params.fov/params.res),params.R,params.PF,BW);
    
    % >>>>>>>>>> calculate diffusion-encoding duration <<<<<<<<<<
            
    C = [2/3*gamma^2*params.G_max^2,gamma^2*params.G_max^2*(params.RF_refoc_dur*1e-3+dur_to_TE*1e-3),0,-params.b_value];
    root = roots(C);
    params.tdiff    = root(find(real(root)>0))*1000;
    pre_readout     = params.tdiff + params.RF_refoc_dur/2+ceil(params.RF_exc_dur/2);
    pre_readout_new = ceil(pre_readout)-(dur_to_TE-floor(dur_to_TE));
    params.TE       = 2*(dur_to_TE + pre_readout_new);
    
    acq_start_time  = params.TE-dur_to_TE;
    
    params.label=strcat('EPI');
    if(params.PF~=1)
        params.label=strcat('PF-EPI');
    end
    params.readout_dur = params.time_acquisition(end);
    params.time_acquisition = params.time_acquisition+acq_start_time;
end


