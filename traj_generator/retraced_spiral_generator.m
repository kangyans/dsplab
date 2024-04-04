function [k,w,Ptn] = retraced_spiral_generator(Nint,res,fov,RR,traj_type)
%
% Function
% ====================
%   Implement of spiral out/ partially spiral in-out/ fully
%   spiral in-out
%
% Requirement:
% ====================
%   Brian H. vdspiral toolbox required.
%
% Parameters:
% ====================
%   - Nint: number of interleaves
%   - fov : field of view in unit of cm
%   - res : resolution in unit of mm
%   - RR  : factor on how many points = (R*Nx) you want to retraced.
%
% (c) Kang Yan (2023), Univeristy of Virginia.
%

  
% ----------------------------------------------------------------------
% Set parameters for sprial design.
% ----------------------------------------------------------------------
smax   = 15000;
gmax   = 4;
T      = .000004;

res    = res*1e3; % mm
rmax   = 5/res;
fov    = fov*100; % cm
N      = fov*10/res;
imsize = [N,N];
fov1   = [fov -fov/2];

[traj,~,~,~,~,~] = vds(smax,gmax,T,Nint,fov1,rmax);

% ----------------------------------------------------------------------
% spiral type pickup
% ----------------------------------------------------------------------

Nx  = length(traj);
ko  = zeros(length(traj),Nint);
wo  = zeros(size(ko));

for i = 1:Nint
    ktmp(:,i) = traj.*exp(-1i*2*pi*(i-1)/Nint);
end

ko  = ktmp./(2*max(abs(ktmp(:)))); % Normalization to [-.5,.5]
wo  = caldcf(ko,imsize);


Ptn = floor(RR*Nx);
if Ptn == 0 && traj_type == "sp"
    k = ko;
    w = wo;
elseif Ptn == 0 && traj_type == "spoi"
    % spiral out-in
    k  = [ko; flip(ko(2:end,:))];
    w  = [wo; flip(wo(2:end,:))];

    % spiral in-out
    % k  = [flip(ko(2:end,:)); ko];
    % w  = [flip(wo(2:end,:)); wo];

elseif Ptn > 0
    kr = ko(1:Ptn,:);
    wr = wo(1:Ptn,:);
    k  = [flip(kr); ko(2:end,:); flip(ko(Ptn+1:end,:))];
    w  = [flip(wr); wo(2:end,:); flip(wo(Ptn+1:end,:))];
end

end
