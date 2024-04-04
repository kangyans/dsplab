function fmap = off_res_fmap_sim(imsize, fscale)
sigma  = imsize(1)/2;
[X, Y] = meshgrid(1:imsize(1), 1:imsize(2));
fmap   = (1-exp(-1/(sigma^2)*((Y-imsize(1)/2).^2 + (X-imsize(2)/2).^2)))*fscale;
end