% ----------------------------------------------------------------------- %
%  
% Fourth order diffusion problem by Marco Pingaro
%
% solve a(u,v) = F(v)    for all v H_{0}^{2}
%
% the bilinear form is a(u,v) = A(x)*D^2 u : D^2 u + B(x)*grad u * grad v
% Omogeneus Dirichlet condition.
%
clear all; close all; clc;
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
base = 1; altezza = 1;
F = 1;                           % Carico
A = 10.0;
B = 1.0;
% --------------------------------------------------
p11 =[0 0]; p12 =[base 0]; p22 =[base altezza]; p21 =[0 altezza];
srf_i = nrb4surf(p11,p12,p21,p22);
% Raffinamento (k-raffinament)

new_knots = linspace (0, 1, 21);
new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf_i, {new_knots, new_knots});

% Save Geometry
nrbexport(srf_i,'rectangular.txt');
problem_data.geo_name = 'rectangular.txt';
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  =[];          % Define Neumann conditions
problem_data.drchlt_sides_u=[1 3];     % Define Dirichlet conditions
problem_data.drchlt_sides_r=[];        % Bordi incastrati
% Physical parameters
problem_data.c_diff  = @(x, y) A*ones(size(x));        % Laplace
problem_data.d_diff  = @(x, y) B*ones(size(x));        % Bilaplace
% Source and boundary terms
problem_data.f = @(x, y) F*ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];  % Degree of the splines
method_data.regularity = [2 2];  % Regularity of the splines
method_data.nsub       = [21 21];  % Number of subdivisions
method_data.nquad      = [4 4];  % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] = ...
              solve_fourth_order_diffusion_iso (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FourthOrderDiffusion';
vtk_pts = {linspace(0, 1, 41), linspace(0, 1, 41)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% PLOT IN MATLAB
% Plottagio geometria e control points
figure
nrbctrlplot(srf_i)
figure
nrbkntplot(srf_i)
% Plottaggio raffinamento
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
figure
contour (X, Y, eu)