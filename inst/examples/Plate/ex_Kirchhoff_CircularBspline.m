%--------------------------------------------------------------------------
%-------- EXERCISE: KIRCHHOFF CIRCULAR PLATE (B-spline)
%--------------------------------------------------------------------------
clear all; close all; clc;
%-------------------------------------------------
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
Radius = 1.0;                     % Radius of Plate
% Costanti elastiche
E  =  1;                          % Young modulus   
nu = 0;                           % Poisson modulus
tr = 1;                           % Spessore piastra         
D  = E*tr*tr*tr/(12*(1-nu*nu));   % Flexural rigidity of the plate
p  = -1;                          % Carico distribuito
%--------------------------------------------------
Radius = 1;

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
srf = nrbmak(coefs, knots);
srf = nrbtform(srf, vecscale([Radius Radius 0]));
%
% Raffinamento (k-raffinament)
new_knots = linspace (0, 1, 11);
new_knots = new_knots (2:end-1);
srf_refined = nrbkntins(srf, {new_knots, new_knots});
% Save Geometry
nrbexport(srf,'plate_Circular_KirchhoffBspline.txt');
problem_data.geo_name = 'plate_Circular_KirchhoffBspline.txt';
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  =[];            % Define Neumann conditions
problem_data.drchlt_sides_u=[1 2 3 4];   % Define Dirichlet conditions
problem_data.drchlt_sides_r=[1 2 3 4];   % Bordi incastrati
% Physical parameters
problem_data.c_diff  = @(x, y) D*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) p*ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];  % Degree of the splines
method_data.regularity = [2 2];  % Regularity of the splines
method_data.nsub       = [41 41];% Number of subdivisions
method_data.nquad      = [4 4];  % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_2d_iso (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN
%
% EXPORT TO PARAVIEW
output_file = 'PlateCircularKirchhoffNURBS';
vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
%
% PLOT IN MATLAB
% Plottagio geometria e control points
figure
nrbctrlplot(srf)
figure
nrbkntplot(srf)
% Plottaggio raffinamento
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
% Max Deflection
max_spost = min(min(eu))
soluz_anal = 3*(5+nu)*(1-nu)*p*Radius*Radius*Radius*Radius/(16*E*tr*tr*tr)
%-----------------------------------------------------