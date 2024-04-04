% A simple way to estimate b-value with given gradient.

clear
clc
close all


gamma = 267.522e6;

g_max = 30;


dt1 = 24;
dt2 = 12;
dt3 = 3;
dt4 = 12;
dt5 = 0;

t1 =  1:dt1;
t2 = dt1:dt1+dt2;
t3 = dt1+dt2:dt1+dt2+dt3;
t4 = dt1+dt2+dt3:dt1+dt2+dt3+dt4;
t5 = dt1+dt2+dt3+dt4:dt1+dt2+dt3+dt4+dt5;

t = [0 t1 t2 t3 t4 t5]*1e-3;

g1 =  g_max *ones(1,length(t1));
g2 = -g_max *ones(1,length(t2));
g3 = -    0 *ones(1,length(t3));
g4 = -g_max *ones(1,length(t4));
g5 =      0 *ones(1,length(t5));

g = [g1(1) g1 g2 g3 g4 g5]*1e-3;

k      =  gamma*cumsum(g)/1e6;
btmp   =  cumsum(k.^2)/1e3;
b      =  btmp(end)


figure(1),
subplot(211),plot(t,g);

ylabel('Diffusion Gradient');
xlabel('Time / s');

subplot(212),plot(t, k, 'r--');

ylabel('K-space');
xlabel('Time / s');