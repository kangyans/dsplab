function I = csphantom (n)
%%CSPHANTOM  Test image for compressive sensing reconstructions.
%
%  IMG = CSPHANTOM(N) returns an NxN phantom with various features designed
%  to fail a TV-regularized CS MRI reconstruction.  
%
%  Written by David S. Smith, 30 Oct 2010.
%
%  Version 1.0

if n < 8, error('Smallest sensible phantom is 8 x 8.'); end
n0 = n;
n = n - round(0.1*n);
bg = 0.5;

% quadrant size
nq = round(n/2);
[X,Y] = meshgrid(1:nq, 1:nq);

% QUADRANT II: low-contrast circles
nc = 4;  % # rows and columns of circles
Q1 = bg*ones(nq,nq);
inc = nq / (nc+1);
for j = 1:nc
  rad = nq / (11 + j^2);
  for k = 1:nc
    val  = bg - 0.1/k;
    rctr = round(j*inc);
    cctr = round(k*inc);
    Z2 = (X-cctr).^2 + (Y-rctr).^2 <= rad^2;
    Q1(Z2) = val;
  end
end

% QUADRANT I: diagonal ramp
Q2 = bg*ones(nq,nq);
scale = nq/4;
ctr = round(nq/2);
Z = abs(X-ctr) + abs(Y-ctr);
Z1 = Z <= scale;
Z2 = Z > scale & Z <= 2*scale;
Q2(Z1) = Q2(Z1) * 2 .* (1 - Z(Z1) / scale);
Q2(Z2) = Q2(Z2) .* (Z(Z2) / scale - 1);
d = round(3*nq/16);
r = [d nq-d];
c = [d nq-d];
w = round(d/4);
for j = 1:2
  for k = 1:2
   Q2(r(j)-w:r(j)+w,c(k)-w:c(k)+w) = bg*(1+j/37+k/43);
  end
end


% QUADRANT III: Gaussian bumps and quadratic hole
ctr = round(nq/2);
scale = nq / 4;
R2 = (X-ctr).^2 + (Y-ctr).^2;
%Q3 = bg * (1 - exp(- R2 / scale^2));
% Q3(R2 <= scale^2) = bg * R2(R2 <= scale^2) / scale^2;
Q3 = bg * R2 / scale^2;
Q3(Q3 > bg) = bg;
ctr = round([nq/5 4*nq/5]);
scale = scale / 6;
for r = 1:2
  for c = 1:2
    Q3 = Q3 + bg * exp(-((X-ctr(r)).^2+(Y-ctr(c)).^2)/scale^2);
  end
end

% QUADRANT IV: line pairs
Q4 = bg*ones(nq,nq);
offset = round(nq/8);
w = max(1,offset);
c = max(1,offset);
r = max(1,offset);
d = 2;
while r <= nq - offset  % horizontal lines
  Q4(r,c:c+w) = 0;
  r = r + d;
  d = d + 1;
end
d = 2;
c = c + offset + w;
r = max(1,offset);
while c <= nq - offset  % vertical lines
  Q4(r:r+w,c) = 0;
  c = c + d;
  d = d + 1;
end
ctr = round(nq/2) + offset;               % concentric circles
R = sqrt((X-ctr).^2 + (Y-ctr).^2);
rad = 2;
thick = 0.5;
while rad < (nq/2 - offset)
  rim = abs(R - rad) <= thick;
  Q4(rim) = 0;
  rad = round(1.5*rad + 1.5);
end

npadf = (n0 - 2*nq) / 2;
npad = round(npadf);

I = padarray([Q1 Q2; Q3 Q4], [npad npad]);
if npadf ~= npad
  I = I(1:end-1,1:end-1); % cut off last column and row if odd dimensions
end

