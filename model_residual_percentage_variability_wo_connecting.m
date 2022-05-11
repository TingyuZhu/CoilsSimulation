%% residue of the model
clearvars
close all
clc

% side wall
% 4 turns

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

r = 1

x = r*linspace(-0.5,0.5,20);    % x [m]
y = r*linspace(-0.5,0.5,20);    % y [m]
z = r*linspace(-0.5,0.5,20);    % z [m]

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

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

figure()
plot(model1)
title 'side wall model&B'
hold on
plot(B)

figure()
plot(residual1)
title 'side wall residual'

100*std(residual1)/std(B) % percentage of variability not fit the real data
100*max(abs(residual1))/max(abs(B))


% second order
C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

figure() 
plot(model2) 
hold on 
plot(B)

figure()
plot(residual1)

100*std(residual2)/std(B) % percentage of variability not fit the real data
100*max(abs(residual2))/max(abs(B))

% third order
C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

figure()
plot(model3)
hold on
plot(B)

figure()
plot(residual3)


100*std(residual3)/std(B) % percentage of variability not fit the real data
100*max(abs(residual3))/max(abs(B))





%% residue of the model
clearvars
close all
clc

% front/back
% 6 turns

BSmag = BSmag_init(); % Initialize BSmag analysis

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

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

figure()
plot(model1)
title 'side wall model&B'
hold on
plot(B)

figure()
plot(residual1)
title 'side wall residual'

% second order
C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

figure() 
plot(model2) 
hold on 
plot(B)

figure()
plot(residual1)

% third order
C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

figure()
plot(model3)
hold on
plot(B)

figure()
plot(residual3)

std(residual3)/std(model3) % variability not fit the model


%% residue of the model
clearvars
close all
clc

% ceiling/floor
% 2 turns

BSmag = BSmag_init(); % Initialize BSmag analysis

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

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

figure()
plot(model1)
title 'side wall model&B'
hold on
plot(B)

figure()
plot(residual1)
title 'side wall residual'

% second order
C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

figure() 
plot(model2) 
hold on 
plot(B)

figure()
plot(residual1)

% third order
C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

figure()
plot(model3)
hold on
plot(B)

figure()
plot(residual3)

std(residual3)/std(model3) % variability not fit the model


%% residue of the model
clearvars
close all
clc

% cube - side
% 2 turns

BSmag = BSmag_init(); % Initialize BSmag analysis

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


norms = zeros(10,3);

r = 0.4
x = r*linspace(-0.5,0.5,20);    % x [m]
y = r*linspace(-0.5,0.5,20);    % y [m]
z = r*linspace(-0.5,0.5,20);    % z [m]

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

% first order
C1 = pinv(H1)*B;
model1 = H1*C1;
residual1 = B - model1;

figure()
plot(model1)
title 'side wall model&B'
hold on
plot(B)

figure()
plot(residual1)
title 'side wall residual'

% second order
C2 = pinv(H2)*B;
model2 = H2*C2;
residual2 = B - model2;

figure() 
plot(model2) 
hold on 
plot(B)

figure()
plot(residual2)

100*std(residual2)/std(B) % percentage of variability not fit the real data
100*max(abs(residual2))/max(abs(B))

% third order
C3 = pinv(H3)*B;
model3 = H3*C3;
residual3 = B - model3;

figure()
plot(model3)
hold on
plot(B)

figure()
plot(residual3)

100*std(residual3)/std(B) % percentage of variability not fit the real data
100*max(abs(residual3))/max(abs(B)) % max of error for the real data

