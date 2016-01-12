%--------------------------------------------------------------------------
%-------- EXERCISE: KIRCHHOFF PLATE (NURBS)
%
%
%--------------------------------------------------------------------------
clear; 
close all; 
clc;
%--------------------------------------------------------------------------
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
base = 1; altezza = 1;
% Costanti elastiche
E  =  1;                           % Young modulus   
nu = 0.0;                          % Poisson modulus
%tr = 0.2;                         % Spessore piastra         
%D  = E*tr*tr*tr/(12*(1-nu*nu));   % Flexural rigidity of the plate
D = 1;
p = 0;                             % carico distribuito
%--------------------------------------------------
p11 =[0 0]; p12 =[base 0]; p22 =[base altezza]; p21 =[0 altezza];
srf_i = nrb4surf(p11,p12,p21,p22);
% Raffinamento (k-raffinament)
new_knots = linspace (0, 1, 3);
new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf_i, {new_knots, new_knots});

% Save Geometry
nrbexport(srf_i,'plate_KirchhoffClamped.txt');
problem_data.geo_name = 'plate_KirchhoffClamped.txt';
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides  = [];          % Define Neumann conditions
problem_data.drchlt_sides_u= [1 2];     % Define Dirichlet conditions u
problem_data.drchlt_sides_r= [1 2] ;    % Define Dirichlet condition du/dn

% Physical parameters
problem_data.c_diff  = @(x, y) D*ones(size(x));
problem_data.d_diff  = @(x, y) D*nu*ones(size(x));
problem_data.e_diff  = @(x, y) 2*D*(1-nu)*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) p*ones(size(x));
%problem_data.h = @boundary_plate_h_drchlt;
%problem_data.g = @(x, y, ind) -ones(size(x));

problem_data.h = @(x, y, ind) zeros(size(x));
problem_data.g = @boundary_plate_g_drchlt;
%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [2 2];   % Degree of the splines
method_data.regularity = [1 1];   % Regularity of the splines
method_data.nsub       = [10 10]; % Number of subdivisions
method_data.nquad      = [3 3];   % Points for the Gaussian quadrature rule
%-------------------------------------------------
% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_plate_kirchhoff (problem_data, method_data);

    % solve_bilaplace_2d_NURBS_iso (problem_data, method_data);
    % solve_bilaplace_GRADGRAD_2d_iso (problem_data, method_data);
    % solve_plate_kirchhoff (problem_data, method_data);
%-------------------------------------------------
% POST-PROCESSIN
disp(msh.nel);  % numero elementi usato!
% EXPORT TO PARAVIEW
output_file = 'PlateKirchhoffClamped';
vtk_pts = {linspace(0, 1, 101), linspace(0, 1, 101)};
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
title ('Numerical solution','fontsize',30)
axis equal


nurbs = geometry.nurbs;
weights = nurbs.coefs(4,:,:);

nurbs.coefs(3,:,:) = nurbs.coefs(3,:,:) + reshape(u,size(weights)) .* weights;
nrbctrlplot(nurbs)
% def_geom = geo_deform (u, space, geometry);
% nrbplot (def_geom.nurbs, [20 20], 'light', 'on')
% view(2)
% title ('Deformed configuration')

% Max Deflection
format long
max_spost  = max(max(eu));
soluz_anal = 0.00126*p*min(base,altezza)^4/D; % clamped
%soluz_anal = 0.004062352626538; % simply supported

fprintf('Soluzione computata = %s \n',max_spost);
fprintf('Soluzione analitica = %s \n',soluz_anal);
%-----------------------------------------------------

%% Soluzione analitica:
x = 0:0.001:base;
sol = (pi/4/base).*x.^3-(pi/4).*x.^2;

le(:,1) = X(:,101);
le(:,2) = eu(:,101);

figure, 
plot(le(:,1), le(:,2),'-o')
hold on
plot(x, sol)
hold off


