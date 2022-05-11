clearvars
close all
clc

% make a rectangular coil with size of 3.2 m by 2 m.
% cube

BSmag = BSmag_init(); % Initialize BSmag analysis

% Define the corners of the loop (including returning to the first coil!)
% Note, units are [m]
% dimension [x y z] = [3.2 4.2 2.0]

%% cube side walls

% anti-clockwise, 2-turns
Gamma7 = [-1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
          -1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
          -1, -1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma7,I,dGamma);

Gamma8 = [1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
          1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
          1, -1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma8,I,dGamma);

% Generate the field points (where we want to see field values)
% centre meter cube

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z); 

% BSmag_plot_field_points(BSmag,xM,yM,zM); % -> shows the field point line
sp = spm_mesh_sphere(5); 
v = sp.vertices*.5;

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,v(:,1),v(:,2),v(:,3));      

% Plot B/|B|
figure(1)
    quiver3(X,Y,Z,BX,BY,BZ,'b')

% normalisation
normB=sqrt(BX.^2+BY.^2+BZ.^2);


% norm
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= normB;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

 
 % BX 
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BX;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BY
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BY;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BZ
 f = figure()
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



    
%% cube front door and back wall

clearvars
close all
clc

BSmag = BSmag_init(); % Initialize BSmag analysis

% anti-clockwise, 2-turns
Gamma = [1, -1, -1; 1, -1, 1; -1, -1, 1; -1, -1, -1; 
    1, -1, -1; 1, -1, 1; -1, -1, 1; -1, -1, -1; 
    1, -1, -1; 
    
    1, 1, -1; 1, 1, 1; -1, 1, 1; -1, 1, -1; 
    1, 1, -1; 1, 1, 1; -1, 1, 1; -1, 1, -1; 
    1, -1, -1;] - [0, 0, 0]; 

I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

% Generate the field points (where we want to see field values)
% centre meter cube

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z); 

% BSmag_plot_field_points(BSmag,xM,yM,zM); % -> shows the field point line
sp = spm_mesh_sphere(5); 
v = sp.vertices*.5;

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,v(:,1),v(:,2),v(:,3));      

% Plot B/|B|
figure(1)
    quiver3(X,Y,Z,BX,BY,BZ,'b')

% normalisation
normB=sqrt(BX.^2+BY.^2+BZ.^2);


 % norm
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= normB;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

 
 % BX 
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BX;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BY
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BY;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BZ
 f = figure()
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

%% cube ceiling and floor

clearvars
close all
clc

BSmag = BSmag_init(); % Initialize BSmag analysis

% anti-clockwise, 2-turns
Gamma = [-1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
    -1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
    -1, -1, -1; 
    
    1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
    1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
    1, -1, -1;] - [0, 0, 0]; 

I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

% Generate the field points (where we want to see field values)
% centre meter cube

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z); 

% BSmag_plot_field_points(BSmag,xM,yM,zM); % -> shows the field point line
sp = spm_mesh_sphere(5); 
v = sp.vertices*.5;

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,v(:,1),v(:,2),v(:,3));      

% Plot B/|B|
figure(1)
    quiver3(X,Y,Z,BX,BY,BZ,'b')

% normalisation
normB=sqrt(BX.^2+BY.^2+BZ.^2);


 % norm
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= normB;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])

 
 % BX 
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BX;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BY
 f = figure()
 p= [];
 p.vertices=sp.vertices;
 p.faces=sp.faces;
 p.EdgeColor='none';
 col= BY;
 patch(p,'FaceVertexCData',col,'FaceColor','interp');
 view([59,28])


 % BZ
 f = figure()
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

