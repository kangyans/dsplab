function params = seq_time_design(params)
% ======================================================================
% Define relaxation times 
% ======================================================================
%
% T1, T2, T2* in (ms) used from:
%
%   [1].Cox and Gowland., 2010
%   [2].Peters    et al., 2007
%   [3].Rooney    et al., 2007
%   [4].Wansapura et al., 1999
% ----------------------------------------------------------------------

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

% ======================================================================
% Set relaxation times based on tissue and field strength 
% ======================================================================

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

% ======================================================================
% Trajectory Generation
% ======================================================================

if(params.traj_type == "sp" || params.traj_type == "spoi")
        
    [params.k,params.w]  = retraced_spiral_generator(params.Nint,params.res,params.fov,0,params.traj_type); 
    params.tdiff = diff_design(params.b_value/1e6, params.G_max*1e3,params.dRF_refoc,'kbct','rectangle'); % ms
    params.dt1   = (params.TE/2 - params.tdiff)/2 - (params.dRF_exc/2 - params.dRF_refoc/2);
    params.dt2   = (params.TE/2 - params.tdiff)/2 + (params.dRF_exc/2 - params.dRF_refoc/2);
    minTE  = 2*(params.tdiff + params.dRF_refoc/2 + params.dRF_exc/2);  
    if (params.TE  < minTE)
        error('Please increase TE!!!');
    end
    dwelltime   = 5e-3; % ms, 1 us 
    params.tADC = params.TE + (0:size(params.k,1) - 1)'*dwelltime;   
    params.dRO  = size(params.k,1)*dwelltime;
    params.label   = strcat('Spiral');
    params.kloc_PF = [];

end


if(params.traj_type == "spaio")

    [params.k,params.w,params.extraPtn]  = retraced_spiral_generator(params.Nint,params.res,params.fov,params.extraR,params.traj_type);
    params.tdiff = diff_design(params.b_value/1e6, params.G_max*1e3,params.dRF_refoc,'kbct','rectangle');
    dwelltime    = 5e-3;
    params.dt1   = (params.TE/2 - params.tdiff)/2 - (params.dRF_exc/2 - params.dRF_refoc/2);
    params.dt2   = (params.TE/2 - params.tdiff)/2 + (params.dRF_exc/2 - params.dRF_refoc/2);
    if (params.dt1  < params.extraPtn * dwelltime)
        error('Please increase TE!!!');
    end

    params.tADC = params.TE + (-params.extraPtn:(size(params.k,1)-params.extraPtn) - 1)'*dwelltime;    
    params.readout_dur = size(params.k,1)*dwelltime;
    params.label       = strcat('Spiral Asymmetric in-out');
    params.kloc_PF     = [];
end

end

