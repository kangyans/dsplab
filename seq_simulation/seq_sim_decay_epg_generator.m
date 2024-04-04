function sig_decay = seq_sim_decay_epg_generator(seq)

% **************************************************
% Set parameter
% **************************************************
N     = seq.num_dummy           ;
G     = seq.G_max*10000         ; % G/m
tG    = seq.tdiff*1e-3          ; % s
dt1   = seq.dt1*1e-3            ; % s
dt2   = seq.dt2*1e-3            ; % s 
T1    = seq.T1 * 1e-3           ; % s
T2    = seq.T2 * 1e-3           ; % s
T2s   = seq.T2s* 1e-3           ; % s
TE    = seq.TE * 1e-3           ; % s
TR    = seq.TR * 1e-3           ; % s

gamma = 4258*2*pi               ; % gamma over 2*pi
D     = 1*10^(-9)               ; % Diffusion (m^2/s)
dkn   = gamma*G*tG              ;

% **************************************************
% EPG simulation of Diffusion spin-echo sequence
% **************************************************
    
F = [0;0;1]          ;
F = epg_rf(F,pi/2,0) ;			% 90 RF Rotation.
F = epg_grelax(F,T1,T2s,dt1,dkn,D,0,0);	% Relaxation, Diffusion, "dt1"
F = epg_grelax(F,T1,T2s, tG,dkn,D,1,0);	% Relaxation, Diffusion, Gradient
F = epg_grelax(F,T1,T2s,dt2,dkn,D,0,0);	% Relaxation, Diffusion, "dt2"

F = epg_rf(F,pi,pi/2);          % 180 RF Rotation.
F = epg_grelax(F,T1,T2,dt2,dkn,D,0,0);	% Relaxation, Diffusion, "dt2"
F = epg_grelax(F,T1,T2,tG ,dkn,D,1,0);	% Relaxation, Diffusion, Gradient

time_pre_TE   = seq.tADC(seq.tADC <  seq.TE);
time_post_TE  = seq.tADC(seq.tADC >= seq.TE);

for i = 1: length(time_pre_TE)
    dt = (seq.dt1 - (seq.TE - time_pre_TE))*1e-3;
    F1 = epg_grelax(F,T1,T2,dt(i) ,dkn,D,0,0);
    decay_pre_TE(i) = abs(F1(1,1));
end

F = epg_grelax(F,T1,T2,dt1,dkn,D,0,0);	% Relaxation, Diffusion, "dt1"
for i = 1: length(time_post_TE)
    dt = ( time_post_TE - seq.TE)*1e-3;
    F2 = epg_grelax(F,T1,T2s,dt(i) ,dkn,D,0,0);
    decay_post_TE(i) = abs(F2(1,1));
end

if isempty(time_pre_TE) == 1
    sig_decay = [decay_post_TE(:)];
else
    sig_decay = [decay_pre_TE(:);decay_post_TE(:)];
end

end
