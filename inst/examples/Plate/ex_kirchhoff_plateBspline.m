%--------------------------------------------------------------------------
%-------- EXERCISE: KIRCHHOFF PLATE (B-SPLINE)
%--------------------------------------------------------------------------
clear all; close all; clc;
%-------------------------------------------------
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
base = 1; altezza = 1;
% Costanti elastiche
E  =  30;                         % Young modulus   
nu = 0.3;                         % Poisson modulus
tr = 0.2;                         % Spessore piastra         
D  =   E*tr*tr*tr/(12*(1-nu*nu)); % Flexural rigidity of the plate             
p = -1;                           % Carico distribuito
%--------------------------------------------------
p11 =[0 0]; p12 =[base 0]; p22 =[base altezza]; p21 =[0 altezza];
srf = nrb4surf(p11,p12,p21,p22);
% Raffinamento (k-raffinament)
new_knots = linspace (0, 1, 11);
new_knots = new_knots (2:end-1);
srf_refined = nrbkntins(srf, {new_knots, new_knots});
% Save Geometry
nrbexport(srf,'plate_Kirchhoff.txt');
problem_data.geo_name = 'plate_Kirchhoff.txt';
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  =[];            % Define Neumann conditions
problem_data.drchlt_sides_u=[1 2 3 4];   % Define Dirichlet conditions
problem_data.drchlt_sides_r=[1 2 3 4];   % Bordi incastrati 
% Physical parameters
problem_data.c_diff  = @(x, y) D*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) -ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];  % Degree of the splines
method_data.regularity = [2 2];  % Regularity of the splines
method_data.nsub       = [11 11];% Number of subdivisions
method_data.nquad      = [4 4];  % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_2d_iso (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'PlateKirchhoff';
vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

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
max_spost  = min(min(eu))
soluz_anal = 0.00406*p*min(base,altezza)^4/D
%-----------------------------------------------------