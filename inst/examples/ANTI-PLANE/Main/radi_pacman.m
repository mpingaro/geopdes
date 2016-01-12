% RADI version
% Using equation (5)
clear; 
close all; 
clc

% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KIII = 0.5;   % Se varia cambiare : boundary_pacman_mixed_bc_h_drchlt  !!!!
mu = 1;
eta = 0.0;
l = 0.05;
%Radius = 100*l; % Radius
Radius = 10;
lt = l*sqrt(1+eta);
ls = l/sqrt(2);

c1 = 1.0;
c2 = lt*lt/2;
c3 = ls*ls - lt*lt/2;
points = 500;
% Definizione della geometria
knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
coefs = zeros(4,3,3);
coefs(:, 1, 1) = [1, 0, 0, 1]; 
coefs(:, 2, 1) = [0, 0, 0, 1]; %  0.5 prima
coefs(:, 3, 1) = [0, 0, 0, 1];
coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 2, 2) = [0, 0, 0, 1];
coefs(:, 3, 2) = [0, 0, 0, 1]; % -0.5 prima
coefs(:, 1, 3) = [0, 1, 0, 1]; 
coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 3) = [-1, 0, 0, 1];

srf = nrbmak(coefs, knots);
srf = nrbtform(srf, vecscale([Radius Radius 0]));

% % PLOT IN MATLAB.
% % Plottagio geometria e control points
% figure
% nrbctrlplot(srf)
% % Plottaggio raffinamento
% figure
% nrbkntplot(srf)

% Elevazione d'ordine
% srf = nrbdegelev(srf,[0 1]);
% k-raffinament
new_knots = linspace (0, 1, 10);
%new_knots = [0,logspace(-3,0,10)];
%new_knots = [ones(1,10)-logspace(0, -3,10),1];

new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf, {new_knots, new_knots});

% Secondo modo---
% [~,~,new_knots] = kntrfine(srf.knots,[9 9],[2 2],[1 1]);
% srf_refined = nrbkntins(srf,new_knots);
% Save Geometry
%nrbexport(srf,'pacman.txt');
%problem_data.geo_name = 'pacman.txt';
problem_data.geo_name = srf_r;
% ----------------------------------------------------------------------- %
% Type of boundary conditions for each side of the domain
% 1 Lato arco destro            |---> u_d = u_far;
% 2 Lato orizzontale sinistro;  |---> u libero;
% 3 Lato orizzontale destro;    |---> u_d = 0;
% 4 Lato arco sinistro          |---> u_d = u_far;

problem_data.nmnn_sides  = [];           % Define Neumann conditions
problem_data.drchlt_sides_u= [1, 3, 4];  % Define Dirichlet conditions u
problem_data.drchlt_sides_r= [];  % Define Dirichlet condition du/dn

% ----------------------------------------------------------------------- %
% Physical parameters

problem_data.c_diff  = @(x, y) c1.*ones(size(x));        % Grad x Grad
problem_data.d_diff  = @(x, y) c2.*ones(size(x));        % GradGrad x GradGrad
problem_data.e_diff  = @(x, y) c3.*ones(size(x));        % Laplace x Laplace

% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));                 % Bulk
problem_data.g = @boundary_pacman_mixed_bc_g_drchlt;     % Rotation
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;     % Displacement

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [srf_r.dim-1 srf_r.dim-1];       % Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];       % Regularity of the splines
method_data.nsub       = [1 1];                           % Number of subdivisions
method_data.nquad      = [3 3];                           % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] = solve_fracture_radi (problem_data, method_data);

%% POST-PROCESSIN
% EXPORT TO PARAVIEW
output_file = 'FractureRadi';
vtk_pts = {linspace(0, 1, points), linspace(0, 1, points)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
% Plottagio geometria e control points
% figure
% nrbctrlplot(srf)
% % Plottaggio raffinamento
% figure
% nrbkntplot(srf)
% 
% figure,
% nrbctrlplot(srf_r)
% 
% figure
% nrbkntplot(srf_r)

% Solution
figure,
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
%leap_SG_COMP(:,1) = [X(:,1);X(points,:)'];
%leap_SG_COMP(:,2) = [eu(:,1);eu(points,:)'];
leap_SG_COMP(:,1) = X(points,:);
leap_SG_COMP(:,2) = eu(points,:);
anal_sol = sqrt(2/pi)*sqrt(-leap_SG_COMP(:,1)).*KIII./mu;
%
figure
plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),leap_SG_COMP(:,1),anal_sol)
axis equal
%print('lip2.eps')
%save leap_SG_COMP
% ----------------------------------------------------------------------- %


% figure, plot(theta,w1*G*sqrt(pi/(2*l))/KII)
% title ('FIGURE 11 a) Normalized angular variation of the out-of-plane displacement field w at distance 0.1l','fontsize',15)
% xlabel('\theta','fontsize',20)
% ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
% axis([0 180 0 0.4])
% grid on
% 
% figure, plot(theta,w2*G*sqrt(pi/(2*l))/KII)
% title ('FIGURE 11 b) Normalized angular variation of the out-of-plane displacement field w at distance 0.4l','fontsize',15)
% xlabel('\theta','fontsize',20)
% ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
% axis([0 180 0 0.8])
% grid on

% Valutazione flesso sulla lip
% p=polyfit(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),9);
% y=polyval(p,leap_SG_COMP(npt:end,1));
% figure, plot(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),leap_SG_COMP(npt:end,1),y);
% 
% dp=polyder(p) ;
% ddp=polyder(dp);
% flex=roots(ddp);
% disp(abs(max(flex)))
% ----------------------------------------------------------------------- %
figure
plot(leap_SG_COMP(:,1)./l,leap_SG_COMP(:,2).*(mu*sqrt(pi/(2*l))/KIII),'- r', 'linewidth',2 );
title ('FIGURE 9 a) Variation of crack face sliding displacement w along the crack face','fontsize',15)
xlabel('x_{1}/l','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([-5 0 0 2.5])
grid on

figure
plot(log(-leap_SG_COMP(:,1)./l), log(leap_SG_COMP(:,2).*(mu*sqrt(pi/(2*l))/KIII)),'- r', 'linewidth',2 );
title ('FIGURE 9 b) Variation of crack face sliding displacement w along the crack face in logarithmic scales','fontsize',15)
xlabel('log_{10}(x_{1}/l)','fontsize',20)
ylabel('log_{10}[wG(\pi/2l)^{1/2}/K_{III}]','fontsize',20)
axis([-3 1 -3 1])
grid on

figure,
plot(leap_SG_COMP(:,1)./l,anal_sol.*(mu*sqrt(pi/(2*l))/KIII))
title ('First gradient normalized','fontsize', 15)