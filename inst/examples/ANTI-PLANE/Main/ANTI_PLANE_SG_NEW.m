clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
%KII = 1;
%mu = 1;
Ls = 0.2; Lt = 0.2;

C1 =Ls*Ls;
C2 =Ls*Ls-Lt*Lt/2;
C3 =(Ls*Ls)/(Lt*Lt);
% Definizione della geometria
Radius = 1;

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
% srf = nrbdegelev(srf,[0 1]);
% k-raffinament
new_knots = linspace (0, 1, 60);
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

problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1,3,4];
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) C1*ones(size(x));     % Laplace
% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [srf_r.dim-1 srf_r.dim-1];      % Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];      % Regularity of the splines
method_data.nsub       = [1 1];                          % Number of subdivisions
method_data.nquad      = [4 4];                          % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_2d_GRAD_FRAC_iso (problem_data, method_data);
% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, 60), linspace(0, 1, 60)};
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
leap_SG(:,1) = [X(:,1);X(60,:)'];
leap_SG(:,2) = [eu(:,1);eu(60,:)'];
figure
plot(leap_SG(:,1),leap_SG(:,2))
save leap_SG
% ----------------------------------------------------------------------- %