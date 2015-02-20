% EX_LAPLACE_ISO_PLATE: solve the Poisson problem in one quarter of a plate with a hole, discretized with NURBS (isoparametric approach).
clear all; close all; clc;
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data 
% Physical domain, defined as NURBS map given in a text file
p11 =[-1 0]; p12 =[1 0]; p22 =[1 1]; p21 =[-1 1];
srf = nrb4surf(p11,p12,p21,p22);

problem_data.geo_name = srf;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.f = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [2 2];       % Degree of the splines
method_data.regularity = [1 1];       % Regularity of the splines
method_data.nsub       = [20 20];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_laplace_2d_iso (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

%output_file = 'Plate_NRB_Deg3_Reg2_Sub8';

vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
%fprintf ('The result is saved in the file %s \n \n', output_file);
%sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution'), axis tight