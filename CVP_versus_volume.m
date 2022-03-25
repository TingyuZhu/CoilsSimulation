%% residue of the model
clearvars
close all
clc

% mean / standard deviation vs volume
% rectangular
% side walls

BSmag = BSmag_init(); % Initialize BSmag analysis
Gamma = [-1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1;  
    -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1;  
    -1.6, -2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1; 
    1.6, -2.1, -1] - [0, 0, 0]; 
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

norms = zeros(10,3);

r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);

BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   

figure(1)
normB=sqrt(BX.^2+BY.^2+BZ.^2);
quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')

% mean of the normB
mean_of_norm = mean(normB(:));

% stand.devi.
standard_deviation = std(normB(:));

B = [std(BX), std(BY), std(BZ)];

norms(i,1) = mean_of_norm/standard_deviation

end


% rectangular
% door and back wall

BSmag = BSmag_init(); % Initialize BSmag analysis

% anti-clockwise, 4-turns
Gamma = [1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
    1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
    1.6, -2.1, -1; 
    
    1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
    1.6, 2.1, -1;] - [0, 0, 0]; 
    % four turns each rectangular solenoid
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);

BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   

figure(3)
normB=sqrt(BX.^2+BY.^2+BZ.^2);
quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')

% mean of the normB
mean_of_norm = mean(normB(:));

% stand.devi.
standard_deviation = std(normB(:));

B = [std(BX), std(BY), std(BZ)];

norms(i,2) = mean_of_norm/standard_deviation

end


% seiling and floor

Gamma = [1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
    1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
    1.6, -2.1, 1; 
    
    1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
    1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
    1.6, -2.1, -1;] - [0, 0, 0]; 
    % two turns each rectangular solenoid
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);

BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   

figure(3)
normB=sqrt(BX.^2+BY.^2+BZ.^2);
quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')

% mean of the normB
mean_of_norm = mean(normB(:));

% stand.devi.
standard_deviation = std(normB(:));

B = [std(BX), std(BY), std(BZ)];

norms(i,3) = mean_of_norm/standard_deviation

end


figure(2)
plot(r,norms)
figure(3)
semilogy(r,norms)
title 'mean divided by standard deviation versus edge length'
xlabel 'edge of the cube'
ylabel 'normB/std(normB)'
grid on



















