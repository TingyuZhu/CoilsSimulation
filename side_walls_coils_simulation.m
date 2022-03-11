clearvars
close all
clc

% make a rectangular coil with size of 4.2 m by 2 m.

BSmag = BSmag_init(); % Initialize BSmag analysis

% Define the corners of the loop (including returning to the first coil!)
% Note, units are [m]
% dimension [x y z] = [3.2 4.2 2.0]
Gamma = [-1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    -1.598, -2.1, -1; -1.598, -2.1, 1; -1.598, 2.1, 1; -1.598, 2.1, -1; 
    -1.596, -2.1, -1; -1.596, -2.1, 1; -1.596, 2.1, 1; -1.596, 2.1, -1; 
    -1.594, -2.1, -1; -1.594, -2.1, 1; -1.594, 2.1, 1; -1.594, 2.1, -1; 
    -1.6, -2.1, -1; 
    
    1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1; 
    1.598, -2.1, -1; 1.598, -2.1, 1; 1.598, 2.1, 1; 1.598, 2.1, -1; 
    1.596, -2.1, -1; 1.596, -2.1, 1; 1.596, 2.1, 1; 1.596, 2.1, -1; 
    1.594, -2.1, -1; 1.594, -2.1, 1; 1.594, 2.1, 1; 1.594, 2.1, -1; 
    1.6, -2.1, -1] - [0, 0, 0]; 
    % opposite direction of current flow for the left and right wall 
    % four turns each rectangular solenoid

I = 0.03; % filament current [A]
dGamma = 1e-3; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

% Generate the field points (where we want to see field values)
% centre meter cube

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z); 

BSmag_plot_field_points(BSmag,xM,yM,zM); % -> shows the field point line

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);      

% Plot B/|B|
figure(1)
    quiver3(X,Y,Z,BX,BY,BZ,'b')

% normalisation
normB=sqrt(BX.^2+BY.^2+BZ.^2);
    
    % mean of the normB
    mean_of_norm = mean(normB(:));

    % stand.devi.
    standard_deviation = std(normB(:));





