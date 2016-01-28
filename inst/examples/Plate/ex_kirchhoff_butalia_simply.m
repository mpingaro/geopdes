% EX Kirchhoff.
% Trapezoidal plate simply supported in two opposite edges 
% and uniformaly distributed load.
clear; 
close all; 
clc;
%--------------------------------------------------------------------------
% DEFINIZIONE GEOMETRIA
% Geometria iniziale
a = 4;
b = 4;
beta = [pi/12, pi/6, pi/4, pi/3, 5*pi/12]; % 15, 30, 45, 60, 75
% Costanti elastiche
E  =  8.736e7;                     % Young modulus   
nu = 0.3;                          % Poisson modulus
tr = 0.08;                         % Spessore piastra         
D  = E*tr*tr*tr/(12*(1-nu*nu));    % Flexural rigidity of the plate
p = 16/D;                          % carico distribuito
pts = 11;
for ang=beta
%--------------------------------------------------
% Trapezoidal mesh
p11 =[0 0]; p12 =[2*a 0]; p22 =[2*a+2*b*sin(ang) 2*b*cos(ang)]; 
p21 =[2*b*sin(ang) 2*b*cos(ang)];
srf_i = nrb4surf(p11,p12,p21,p22);
% Raffinamento (k-raffinament)
new_knots = linspace (0, 1, 3);
new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf_i, {new_knots, new_knots});

% Save Geometry
nrbexport(srf_i,'plate_KirchhoffSimply.txt');
problem_data.geo_name = 'plate_KirchhoffSimply.txt';
%--------------------------------------------------
% BOUNDARY CONDITIONS 
problem_data.nmnn_sides     = [];        % Define Neumann conditions
problem_data.bound_sides    = [1,2,3,4]; % Define free edges (Boundary integrals)
problem_data.drchlt_sides_u = [3, 4];    % Define Dirichlet conditions u
problem_data.drchlt_sides_r = [];        % Define Dirichlet condition du/dn

% Physical parameters
problem_data.poisson = nu;
problem_data.c_diff  = @(x, y) ones(size(x));
problem_data.d_diff  = @(x, y) nu*ones(size(x));
problem_data.e_diff  = @(x, y) 2*(1-nu)*ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) p*ones(size(x));
%problem_data.h = @boundary_plate_h_drchlt;
%problem_data.g = @boundary_plate_g_drchlt;
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

%--------------------------------------------------
% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [2 2];   % Degree of the splines
method_data.regularity = [1 1];   % Regularity of the splines
method_data.nsub       = [8 8];   % Number of subdivisions
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
% EXPORT TO PARAVIEW
output_file = 'PlateKirchhoffClamped';
vtk_pts = {linspace(0, 1, pts), linspace(0, 1, pts)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% PLOT IN MATLAB
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
title ('Deformed configuration')

% Max Deflection
format long
max_spost  = max(max(eu));
fprintf('Soluzione computata = %s \n',max_spost);
x = (0:2*a/(pts-1):2*a)./a;

w_x = eu(:, (pts-1)/2 + 1)./(p*a^4/D); % middle line x-x
w_y = eu((pts-1)/2 + 1, :)./(p*a^4/D); % middle line y-y

% Save results
ang_d= rad2deg(ang);
val = num2str(ang_d);
name = 'trapezoidal_kirchhoff_ang_';

file_name_1 = strcat(name, val, '_middle_x.txt');
f = fopen(file_name_1, 'w');
fprintf(f,'normalized coordinates v.s. normalized transverse displacement \n');
for i=1:pts
    fprintf(f, '%6.4f \t %6.5e \n', x(i)-1, w_x(i));
end
fclose(f);

file_name_2 = strcat(name, val, '_middle_y.txt');
ff = fopen(file_name_2, 'w');
fprintf(ff,'normalized coordinates v.s. normalized transverse displacement \n');
for i=1:pts
    fprintf(ff, '%6.4f \t %6.5e \n', x(i)-1, w_y(i));
end
fclose(ff);

end