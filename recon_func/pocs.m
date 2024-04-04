
function im = pocs(seq,Niter)
% Parameter:
%   - ksp  : k-space data
%   - Niter: number of iterations
%   - pf   : extra ratio of sampling 
%   - traj : trajectory
%   - dcf  : density compensation function
%
% TODO:
%   - Maybe Hanning windown is required to avoid Gribbs ring | phase.
%   - Make sure your POCS is RIGHT & Good enough!!!!!
%   - POCS + PI ?
%   - Look for Berkin's codes on PI for more reference 
  
% get center data
np      = seq.extraPtn;
[nx,ny] = size(seq.kdata);
nx      = nx -np + 1;
ksp     = zeros(2*nx-1,ny);
k_low   = ksp;

ksp(nx-np+1:end,:) = seq.kdata;
k_low(nx-np+1:nx+np-1,:) = seq.kdata(1:2*np-1,:);

k = cat(1,flip(seq.k(np+1:end,:)), seq.k(np:end,:));
w = cat(1,flip(seq.w(np+1:end,:)), seq.w(np:end,:));

F = NUFFT(k,1,seq.imsize);

theta = angle(F'*(w.*k_low));
figure(111),imagesc(theta);

for i = 1:Niter
    im = F'*(w.*ksp);
    figure(222),imshow(abs(im),[])
    
    % Updata image with phase calculated from low resolution.
    im = abs(im).*exp(1i*theta);

    % Updata k-space.
    k_new   = F*im;
    ksp(1:nx-np,:) = k_new(1:nx-np,:);
    

    %theta = angle(F'*(w.*ksp));
    % Goto next iterative.
end

im = F'*(w.*ksp);

end