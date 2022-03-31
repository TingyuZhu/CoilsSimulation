%% residue of the model
% cube

clearvars
close all
clc

% mean / standard deviation vs volume
% cube

BSmag = BSmag_init(); % Initialize BSmag analysis
Gamma = [-1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
    -1, -1, -1; -1, -1, 1; -1, 1, 1; -1, 1, -1; 
    -1, -1, -1; 
    1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
    1, -1, -1; 1, -1, 1; 1, 1, 1; 1, 1, -1; 
    1, -1, -1;] - [0, 0, 0]; 
I = 0.03; % filament current [A]
dGamma = 1e-1; % filament max discretization step [m]      
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);

norms = zeros(10,3);

r = 0.1:0.1:1;

V = zeros(10,3);

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

V(:,1) = (2*r(i)*linspace(0.0263,0.5,10)).^3;

end


figure(2)
plot(r,norms)
figure(3)
semilogy(r,norms)
title 'mean divided by standard deviation versus edge length'
xlabel 'edge of the cube'
ylabel 'normB/std(normB)'

% versus volume
figure(4)
plot(V,norms);
title 'Mean of norm/Standard deviation vs. Volume'
xlabel 'volume size'
ylabel 'mean(normB(:))/std(normB(:))';

figure(5)
semilogy(V,norms)
title 'Mean of norm/Standard deviation vs. Volume in log scale'
xlabel 'log volume size'
ylabel 'log mean(normB(:))/std(normB(:))';

grid on


