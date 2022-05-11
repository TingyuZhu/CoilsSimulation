%% first order 
% percentage of max of error with the real data or model

clearvars
close all
clc

% Rectangular - side wall pair 
% 4-turns

BSmag = BSmag_init(); % Initialize BSmag analysis
Gamma1 = [-1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma1,I,dGamma);

Gamma2 = [1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma2,I,dGamma);

variability = zeros(10,4);
error = zeros(10,4);
r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

S=[];
S.li=1;
S.reg=1;
S.v=pos;
S.o=or;
H1= spm_opm_vslm(S);

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

variability(i,1) = 100*std(residual1)/std(B); % percentage of variability not fit the real data
error(i,1) = 100*max(abs(residual1))/max(abs(B));

end


% Rectangular - front and back pair 
% 6-turns

Gamma3 = [1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma3,I,dGamma);

Gamma4 = [1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma4,I,dGamma);

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

S=[];
S.li=1;
S.reg=1;
S.v=pos;
S.o=or;
H1= spm_opm_vslm(S);

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

variability(i,2) = 100*std(residual1)/std(B); % percentage of variability not fit the real data
error(i,2) = 100*max(abs(residual1))/max(abs(B));

end



% Rectangular - floor and ceiling pair 
% 2-turns

Gamma5 = [1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma5,I,dGamma);

Gamma6 = [1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma6,I,dGamma); 

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

S=[];
S.li=1;
S.reg=1;
S.v=pos;
S.o=or;
H1= spm_opm_vslm(S);

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

variability(i,3) = 100*std(residual1)/std(B); % percentage of variability not fit the real data
error(i,3) = 100*max(abs(residual1))/max(abs(B));

end



% Cube - side
% 2-turns

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

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

S=[];
S.li=1;
S.reg=1;
S.v=pos;
S.o=or;
H1= spm_opm_vslm(S);

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

variability(i,4) = 100*std(residual1)/std(B); % percentage of variability not fit the real data
error(i,4) = 100*max(abs(residual1))/max(abs(B));

end




figure()
plot(r,variability);
title 'percentage of variability not fit the real data versus edge length when L = 1'
xlabel 'edge length [m]'
ylabel 'percentage of variability [%]'
legend ('side wall','front and back','ceiling and floor','cube')
figure()
plot(r,error);
title 'percentage of max error versus edge length when L = 1'
xlabel 'edge length [m]'
ylabel 'percentage of max error [%]'
legend ('side wall','front and back','ceiling and floor','cube')




%% second order
% percentage of max of error with the real data or model

clearvars
clc


% Rectangular - side wall pair 

BSmag = BSmag_init(); % Initialize BSmag analysis
Gamma1 = [-1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma1,I,dGamma);

Gamma2 = [1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma2,I,dGamma);

variability = zeros(10,4);
error = zeros(10,4);
r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% second order

S=[]
S.li=2;
S.reg=1;
S.v=pos;
S.o=or;
H2= spm_opm_vslm(S);

C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

variability(i,1) = 100*std(residual2)/std(B); % percentage of variability not fit the real data
error(i,1) = 100*max(abs(residual2))/max(abs(B));

end


% Rectangular - front and back pair 

Gamma3 = [1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma3,I,dGamma);

Gamma4 = [1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma4,I,dGamma);

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% second order

S=[]
S.li=2;
S.reg=1;
S.v=pos;
S.o=or;
H2= spm_opm_vslm(S);

C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

variability(i,2) = 100*std(residual2)/std(B); % percentage of variability not fit the real data
error(i,2) = 100*max(abs(residual2))/max(abs(B));

end



% Rectangular - floor and ceiling pair 

Gamma5 = [1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma5,I,dGamma);

Gamma6 = [1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma6,I,dGamma); 

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% second order

S=[]
S.li=2;
S.reg=1;
S.v=pos;
S.o=or;
H2= spm_opm_vslm(S);

C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

variability(i,3) = 100*std(residual2)/std(B); % percentage of variability not fit the real data
error(i,3) = 100*max(abs(residual2))/max(abs(B));

end



% Cube

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

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% second order

S=[]
S.li=2;
S.reg=1;
S.v=pos;
S.o=or;
H2= spm_opm_vslm(S);

C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

variability(i,4) = 100*std(residual2)/std(B); % percentage of variability not fit the real data
error(i,4) = 100*max(abs(residual2))/max(abs(B));

end


figure()
plot(r,variability);
title 'percentage of variability not fit the real data versus edge length when L = 2'
xlabel 'edge length [m]'
ylabel 'percentage of variability [%]'
legend ('side wall','front and back','ceiling and floor','cube')
figure()
plot(r,error);
title 'percentage of max error versus edge length when L = 2'
xlabel 'edge length [m]'
ylabel 'percentage of max error [%]'
legend ('side wall','front and back','ceiling and floor','cube')




%% third order
% percentage of max of error with the real data or model

clearvars
clc

% Rectangular - side wall pair 

BSmag = BSmag_init(); % Initialize BSmag analysis

Gamma1 = [-1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1; -1.6, -2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          -1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma1,I,dGamma);

Gamma2 = [1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1; 1.6, -2.1, 1; 1.6, 2.1, 1; 1.6, 2.1, -1;
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma2,I,dGamma);

variability = zeros(10,4);
error = zeros(10,4);
r = 0.1:0.1:1;

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% third order

S=[]
S.li=3;
S.reg=1;
S.v=pos;
S.o=or;
H3= spm_opm_vslm(S);

C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

variability(i,1) = 100*std(residual3)/std(B); % percentage of variability not fit the real data
error(i,1) = 100*max(abs(residual3))/max(abs(B));

end


% Rectangular - front and back pair 

Gamma3 = [1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, -2.1, 1; -1.6, -2.1, 1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma3,I,dGamma);

Gamma4 = [1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, 2.1, -1; 
          1.6, 2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma4,I,dGamma);

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% third order

S=[]
S.li=3;
S.reg=1;
S.v=pos;
S.o=or;
H3= spm_opm_vslm(S);

C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

variability(i,2) = 100*std(residual3)/std(B); % percentage of variability not fit the real data
error(i,2) = 100*max(abs(residual3))/max(abs(B));

end



% Rectangular - floor and ceiling pair 

Gamma5 = [1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1; 1.6, 2.1, 1; -1.6, 2.1, 1; -1.6, -2.1, 1; 
          1.6, -2.1, 1;] - [0,0,0];
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma5,I,dGamma);

Gamma6 = [1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1; 1.6, 2.1, -1; -1.6, 2.1, -1; -1.6, -2.1, -1; 
          1.6, -2.1, -1;] - [0,0,0];
[BSmag] = BSmag_add_filament(BSmag,Gamma6,I,dGamma); 

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
BSmag_plot_field_points(BSmag,xM,yM,zM);

[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% third order

S=[]
S.li=3;
S.reg=1;
S.v=pos;
S.o=or;
H3= spm_opm_vslm(S);

C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

variability(i,3) = 100*std(residual3)/std(B); % percentage of variability not fit the real data
error(i,3) = 100*max(abs(residual3))/max(abs(B));

end



% Cube

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

for i = 1:length(r)
x = r(i)*linspace(-0.5,0.5,20);    % x [m]
y = r(i)*linspace(-0.5,0.5,20);    % y [m]
z = r(i)*linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   
B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];

% third order

S=[]
S.li=3;
S.reg=1;
S.v=pos;
S.o=or;
H3= spm_opm_vslm(S);

C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

variability(i,4) = 100*std(residual3)/std(B); % percentage of variability not fit the real data
error(i,4) = 100*max(abs(residual3))/max(abs(B));

end


figure()
plot(r,variability);
title 'percentage of variability not fit the real data versus edge length when L = 3'
xlabel 'edge length [m]'
ylabel 'percentage of variability [%]'
legend ('side wall','front and back','ceiling and floor','cube')
figure()
plot(r,error);
title 'percentage of max error versus edge length when L = 3'
xlabel 'edge length [m]'
ylabel 'percentage of max error [%]'
legend ('side wall','front and back','ceiling and floor','cube')




