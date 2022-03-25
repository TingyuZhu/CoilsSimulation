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

x = linspace(-0.5,0.5,20);    % x [m]
y = linspace(-0.5,0.5,20);    % y [m]
z = linspace(-0.5,0.5,20);    % z [m]

[xM, yM, zM] = meshgrid(x,y,z);
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,xM,yM,zM);   

B = [BX(:);BY(:);BZ(:)];

pos = [X(:),Y(:),Z(:); X(:),Y(:),Z(:); X(:),Y(:),Z(:)];
or =  [repmat([1,0,0],size(BX,1)^3,1);repmat([0,1,0],size(BY,1)^3,1);repmat([0,0,1],size(BZ,1)^3,1)];
S=[]
S.li=1;
S.reg=1;
S.v=pos;
S.o=or;
H1= spm_opm_vslm(S);

S=[]
S.li=2;
S.reg=1;
S.v=pos;
S.o=or;
H2= spm_opm_vslm(S);

S=[]
S.li=3;
S.reg=1;
S.v=pos;
S.o=or;
H3= spm_opm_vslm(S);


C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

figure()
plot(model1)
hold on
plot(B)


figure()
plo(residual1)



C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

figure()
plot(model2)
hold on
plot(B)


figure()
plo(residual1)




C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

figure()
plot(model3)
hold on
plot(B)


figure()
plot(residual3)


std(model3)/std(residual3)