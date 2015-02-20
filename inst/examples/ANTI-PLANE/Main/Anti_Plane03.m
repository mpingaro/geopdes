clear all
close all
clc
% DEFINIZIONE GEOMETRIA (PACMAN)

% Costanti di secondo gradiente
KII = 1;
mu = 1;
Ls = 0.2;
Lt = 0.2;

C1 =Ls*Ls;
C2 =Ls*Ls-Lt*Lt/2;
C3 =(Ls*Ls)/(Lt*Lt);
% Definizione della geometria R=1;
Radius = 3;
crv = nrbline([0 0],[Radius 0]);
srf = nrbrevolve(crv,[0 0 0],[0 0 Radius],pi);
% Elevazione d'ordine
% srf = nrbdegelev(srf,[0 1]);
% k-raffinament
new_knots = linspace (0, 1, 10);
new_knots = new_knots (2:end-1);
srf_refined = nrbkntins(srf, {new_knots, new_knots});

% Secondo modo---
% [~,~,new_knots] = kntrfine(srf.knots,[9 9],[2 2],[1 1]);
% srf_refined = nrbkntins(srf,new_knots);
% Save Geometry
nrbexport(srf,'pacman.txt');

problem_data.geo_name = 'pacman.txt';
% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
% Source and boundary terms (bulk)

problem_data.f = @(x, y) 0*x+0*y-1;
problem_data.h = @(x, y, ind) zeros (size (x));

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [9 9];       % Number of subdivisions
method_data.nquad      = [6 6];       % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace_2d_iso (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
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
nrbctrlplot(srf_refined)
figure
nrbkntplot(srf_refined)
% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
%axis tight