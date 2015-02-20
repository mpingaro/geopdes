%-------------------------------------------------------------------------%
%--------      EXERCISE: CIRCULAR PLATE KIRCHHOFF (NURBS)                 %
%-------------------------------------------------------------------------%
% by Marco Pingaro                                                        %
%                                                                         %
%--------------------------------------------------------------------------
clear all; close all; clc;
%--------------------------------------------------------------------------
% DEFINITION OF GEOMETRY
% Initial geometry
Radius = 1.0;                      % radius of Plate
t = 0.1;                           % thickness
E = 1;                             % young modulus
nu = 0.3;                          % poisson modulus
D = E*t^3/(12*(1-nu^2));           % 
p = D;                             % distributed load
%--------------------------------------------------------------------------
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
% Mesh rafinament (k-raffinament)
%new_knots = linspace (0, 1, 3);
%new_knots = new_knots (2:end-1);
%srf_r = nrbkntins(srf_i, {new_knots, new_knots});

srf_r = srf_i;
problem_data.geo_name = srf_r;
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  = [];            % Define Neumann conditions
problem_data.drchlt_sides_u= [1 2 3 4];   % Define Dirichlet conditions
problem_data.drchlt_sides_r= [1 2 3 4];   % Clamped bounds
% Physical parameters
problem_data.c_diff  = @(x, y) D*ones(size(x));
problem_data.d_diff  = @(x, y) D*nu*ones(size(x));
problem_data.e_diff  = @(x, y) 2*D*(1-nu)*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) p*ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
%method_data.degree     = [srf_r.dim-1 srf_r.dim-1];% Degree of the splines
%method_data.regularity = [srf_r.dim-2 srf_r.dim-2];% Regularity of the splines

method_data.degree     = [3 3];                    % Degree of the splines
method_data.regularity = [2 2];                    % Regularity of the splines
method_data.nsub       = [21 21];                  % Number of subdivisions
method_data.nquad      = [4 4];                    % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_bilaplace_2d_NURBS_iso (problem_data, method_data);

    % solve_bilaplace_2d_NURBS_iso (problem_data, method_data);
    % solve_bilaplace_GRADGRAD_2d_iso (problem_data, method_data);
    % solve_plate_kirchhoff (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN
%
% EXPORT TO PARAVIEW
output_file = 'PlateCircularKirchhoffNURBS';
vtk_pts = {linspace(0, 1, 101), linspace(0, 1, 101)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u');
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
title ('Numerical solution','fontsize',25)
axis equal
%% Max Deflection ---------------------------------------------------------
disp(msh.nel);      % numero elementi usato!
disp(srf_r.dim-1);  % Degree B-Spline
disp(srf_r.dim-2);  % Continuity
max_spost = max(max(eu));
%soluz_anal = p/(64*D)*( (5+nu)/(1+nu) ); % Simply supported
soluz_anal = 1/64; % Clamped
soluz_rap1 = soluz_anal/max_spost;
soluz_rap2 = abs((soluz_anal-max_spost)/max_spost);
fprintf('Computational solution at center of the plate = %s \n', max_spost);
fprintf('Analitic solution at center of the plate = %s \n', soluz_anal);
fprintf('Ratio error = %s \n', soluz_rap1);
fprintf('Percentage error = %s \n', soluz_rap2);
%% ------------------------------------------------------------------------
%.. See the paper "The plate paradox for hard and soft simple support"
%                 I. Babuska, J. Pitkaranta. 1989
%   for major detail.
%%