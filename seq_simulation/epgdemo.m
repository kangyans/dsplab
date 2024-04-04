% Explore the EPG...

clear 
clc
close all

% (1). the relationshio between F_n & M_xy
n       = 1;
spacing = .0001;
z       = 0:spacing:1;
M_xy    = 1 + 1i*0; % along x-axis
F_n     = cumsum(M_xy.*exp(-2*pi*1i*n*z));

x       = real(F_n)*spacing;
y       = imag(F_n)*spacing;
figure('Name','EPG Demo'),plot3(x,y,z)
grid on;
xlabel('Mx','FontWeight','bold');
ylabel('My','FontWeight','bold');
zlabel('Voxel Dimension','FontWeight','bold');

% (2). 
