clearvars
close all
clc

% define dimensional information
% design A: covers all walls with dimensions 4.2m x 3.2m x 2m
x(:,1)=(0:.1:4.2); 
y(1,:)=(0:.1:2); 
z= 0; % first in the z = 0 plane

% define the length of the rectangular loop
% 2ax and 2ay are the dimensions of the rectangular current loop
% one side single loop

ax= 4.2/2; 
ay= 2/2; 

% physical expression of the B_0
% B_0=mu0*I*sqrt(a_x^-2+a_y^-2) [T]

I = 0.1; % filament current [A]
mu0 = 4*pi*1e-7; % vacuum permeability [N/A^2]
B0 = mu0*I*sqrt(ax^-2+ay^-2)/pi;

% For a rectangular current loop that is centered at the origin of
% a Cartesian coordinate system, with coordinates (x, y, z), and that
% resides in the z = 0 plane, the magnetic field components are

Xplus=x+ax;
Xminus=x-ax;
Yplus=y+ay;
Yminus=y-ay;

r1= sqrt(Xplus.^2 + Yplus.^2+z^2);
r2= sqrt(Xminus.^2 + Yplus.^2+z^2);
r3= sqrt(Xplus.^2 + Yminus.^2+z^2);
r4= sqrt(Xminus.^2 + Yminus.^2+z^2);

rho = 1/sqrt(ax^-2+ay^-2);

% a rectangular current loop 
% BX = -rho*B0*z/4[1/r1(r1-y-ay)-1/r2(r2-y-ay)-1/r3(r3-y-ay)+1/r4(r4-y+ay)]
% BY = -rho*B0*z/4[1/r1(r1-x-ax)-1/r2(r2-x-ax)-1/r3(r3-x-ax)+1/r4(r4-x+ax)]
% BZ = -rho*B0*z/4[(x+ax)/r1(r1-y-ay) - (y+ay)/r1(r1-x-ax)
%      -(x-ax)/r2(r2-y-ay) - (y+ay)/r2(r2-x+ax) - (x+ax)/r3(r3-y+ay)
%      -(y-ay)/r3(r3-x-ax) + (a-ax)/r4(r4-y+ay) + (y-ay)/r4(r4-x+ax)]
% magnetic field components are: 





% plot the field points
% []

figure(1), hold on, grid on, box on, axis equal
xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]')
view(3), axis tight





