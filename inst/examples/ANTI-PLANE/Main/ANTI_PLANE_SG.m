clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
% Costanti di secondo gradiente
KII = 0.5;
G   = 1;
l = 0.1;
eta = 0.9;
x1 = 0.1*l;
x2 = 0.4*l;
Radius = 40*l;
npoints = 500;

Ls = l/sqrt(2);
Lt = sqrt(eta+1)*l;
C  = G;
C1 = Ls*Ls;
%C2 = Lt*Lt/2;
C2 = Lt*Lt/2 - C1;

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

% Elevazione d'ordine
srf = nrbdegelev(srf,[0 0]);
% k-raffinament
new_knots = linspace (0, 1, 100);
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

problem_data.nmnn_sides   = 2;  % [2 3]
problem_data.drchlt_sides_u = [1, 3, 4];
problem_data.drchlt_sides_r = [1, 3, 4];

% ----------------------------------------------------------------------- %
% Physical parameters
% Physical parameters
problem_data.c_diff  = @(x, y) C.*ones(size(x)) ;        % Laplace
problem_data.d_diff  = @(x, y) C1.*ones(size(x)) ;       % Bilaplace
problem_data.copp = C2 ;                                 % Bordi
% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;
problem_data.hr =@boundary_pacman_mixed_bc_hr_drchlt;

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [srf_r.dim-1 srf_r.dim-1];      % Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];      % Regularity of the splines
method_data.nsub       = [1 1];                          % Number of subdivisions
method_data.nquad      = [srf_r.dim srf_r.dim];          % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_2d_FRAC_iso (problem_data, method_data);
% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, npoints), linspace(0, 1, npoints)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
% Plottagio geometria e control points
figure
nrbctrlplot(srf)
% Plottaggio raffinamento
figure
nrbkntplot(srf)
figure,
nrbctrlplot(srf_r)
figure
nrbkntplot(srf_r)
% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
%axis tight

% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
leap_SG(:,1) = X(npoints,:);
leap_SG(:,2) = eu(npoints,:);
figure
plot(leap_SG(:,1),leap_SG(:,2))
save leap_SG
% Soluzione 1 Gradiente
anal_sol = (KII/G).*sqrt(-leap_SG(:,1).*(2/pi)) ;
%
% Soluzione 2 Gradiente
CIII = 2.8;
rr = linspace(0,1);
s = CIII.*rr.^(3/2).*( 3./ (16.*(Ls/Lt).^2-3) + 1 ) ;
%
% Lip
d = Lt*sqrt(3)/(16*sqrt(2)) *sqrt((Lt/Ls)^4 -32*(Lt/Ls)^2 +128) ;
disp(d)
% ----------------------------------------------------------------------- %
figure
plot(leap_SG(:,1)./l,leap_SG(:,2).*(G*sqrt(pi/(2*l))/KII),'- r', 'linewidth',2 );
title ('FIGURE 9 a) Variation of crack face sliding displacement w along the crack face','fontsize',15)
xlabel('x_{1}/l','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([-5 0 0 2.5])
grid on

figure
plot(log(-leap_SG(:,1)./l), log(leap_SG(:,2).*(G*sqrt(pi/(2*l))/KII)),'- r', 'linewidth',2 );
title ('FIGURE 9 b) Variation of crack face sliding displacement w along the crack face in logarithmic scales','fontsize',15)
xlabel('log_{10}(x_{1}/l)','fontsize',20)
ylabel('log_{10}[wG(\pi/2l)^{1/2}/K_{III}]','fontsize',20)
axis([-3 1 -3 1])
grid on