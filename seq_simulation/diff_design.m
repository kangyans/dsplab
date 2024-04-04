function y = diff_design(x, G_max,dRF,calc_type,model,dt)
% ********************************************************
% function y = diff_design(x,G_max, dRF, calc_type, model,dt)
%
% (a). Trapezoid model
% --------------------------------------------------------
%                       _______
%    Gss               /       \
%   __________________/         \__________________
%
%             _____   <--G_max      _____
%            /     \               /    |\
%    G_diff /  Gd   \             /   Gd| \
%    ______/         \___________/      |  \_______
%             |<-t-> |<---dRF--->|<-t ->|dt| 
%
% 
% (b). Rectangle model
% --------------------------------------------------------
%
%                     ___________
%    Gss             |           |
%   _________________|           |__________________
%
%          __________ <--G_max    __________
%         |          |           |          |             
%         |  G_diff  |           |  G_diff  |           
%    _____|          |___________|          |_______
%         |<-  t   ->|<---dRF--->|<-   t  ->|
% 
% ********************************************************
%
% Input:
% =================
%   - x         : b-value in s/mm^2 or t_diff/t in ms.
%   - G_max     : maximun diffusion-encoding gradient in mT/m.
%   - dRF       : duration of RF refocusing pulse in ms. [Default: 2.76ms]
%   - calc_type : [1] 'kbct': known b-value calculate t_diff;
%                 [2] 'ktcb': known t_diff calculate b_value;
%   - dt        : ramp up/down time in [ms];
% 
% Output:
% =================
%   - y         : flattop time (t_diff - dt) in ms or b-value in s/mm^2.
%
% (c) Kang Yan, 2023.
%


gamma   = 267.522e6  ; % gyromagnetic ratio in (rad/s/T)

G_max   = G_max*1e-3 ; % T/m
dRF     = dRF*1e-3   ; % s


if strcmp(model,'rectangle') || nargin < 6
    dt   = 0;
else 
    dt   = dt*1e-3; % s
end

switch calc_type

    case 'kbct'

        b_value = x*1e6; % s/m^2
        eqs     = [2/3*gamma^2*G_max^2, gamma^2*G_max^2*(dRF + dt), -gamma^2*G_max^2*dt/6, -b_value + gamma^2*G_max^2*dt^3/30];
        res     = roots(eqs);
        y       = res(real(res) > 0)*1000; % convert to ms;
        y       = y - dt*1000;

    case 'ktcb'
        t_diff = x*1e-3; % s
        y = gamma^2*G_max^2*(t_diff^2*(t_diff+dt+dRF - t_diff/3) - dt^2*t_diff/6 + dt^3/30)/1e6;
end


end