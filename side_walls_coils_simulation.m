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

I = -0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

% Generate the field points (where we want to see field values)
% centre meter cube

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z); 

%BSmag_plot_field_points(BSmag,xM,yM,zM); % -> shows the field point line
sp = spm_mesh_sphere(5);
v = sp.vertices*.5;

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,v(:,1),v(:,2),v(:,3));      

% Plot B/|B|
figure(1)
    quiver3(X,Y,Z,BX,BY,BZ,'b')

% normalisation
normB=sqrt(BX.^2+BY.^2+BZ.^2);


f=figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= normB;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

 
 
f=figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BX;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

 
 f=figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BY;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

  
 f=figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BZ;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])
 
 
    % mean of the normB
    mean_of_norm = mean(normB(:))

    % stand.devi.
    standard_deviation = std(normB(:))





