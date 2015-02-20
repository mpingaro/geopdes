%--------------------------------------------------------------------------
%-------- EXERCISE: KIRCHHOFF CIRCULAR PLATE (GRAD GRAD)
%--------------------------------------------------------------------------
%
%
%
clear all; close all; clc;
%-------------------------------------------------
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
Radius = 1.0;                     % Radius of Plate
D = 1;
p = -1;
%--------------------------------------------------
knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
coefs = zeros(4,3,3);
coefs(:, 1, 1) = [1, 0, 0, 1];
coefs(:, 2, 1) = [1/sqrt(2), -1/sqrt(2), 0 , 1/sqrt(2)];
coefs(:, 3, 1) = [0, -1, 0 ,1 ];
coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 2, 2) = [0, 0, 0, sqrt(2)-1];
coefs(:, 3, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 1, 3) = [0, 1, 0, 1];
coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 3) = [-1, 0, 0, 1];
srf_i = nrbmak(coefs, knots);
srf_i = nrbtform(srf_i, vecscale([Radius Radius 0]));
%
% Raffinamento (k-raffinament)
new_knots = linspace (0, 1, 41);
new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf_i, {new_knots, new_knots});
% Save Geometry

%nrbexport(srf,'plate_Circular_KirchhoffNURBS.txt');
%problem_data.geo_name = 'plate_Circular_KirchhoffNURBS.txt';
problem_data.geo_name = srf_r;
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  =[];          % Define Neumann conditions
problem_data.drchlt_sides_u=[1 2 3 4]; % Define Dirichlet conditions
problem_data.drchlt_sides_r=[]; % Bordi incastrati
% Physical parameters
problem_data.c_diff  = @(x, y) D*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) p*ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [srf_r.dim-1 srf_r.dim-1];% Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];% Regularity of the splines
method_data.nsub       = [1 1];                % Number of subdivisions
method_data.nquad      = [4 4];   % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] = ...
    solve_bilaplace_GRADGRAD_2d_iso (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN
%
% EXPORT TO PARAVIEW
output_file = 'PlateCircularKirchhoffNURBS';
vtk_pts = {linspace(0, 1, 71), linspace(0, 1, 71)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
%
% PLOT IN MATLAB
% Plottagio geometria e control points
figure
nrbctrlplot(srf_i)
figure
nrbkntplot(srf_i)
% Plottaggio raffinamento
figure,
nrbctrlplot(srf_r)
%figure
nrbkntplot(srf_r)
% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
% Max Deflection
max_spost = min(min(eu));
soluz_anal = -5/64;
fprintf('Soluzione computata = %d \n', max_spost);
fprintf('Soluzione analitica = %d \n', soluz_anal);
%-----------------------------------------------------