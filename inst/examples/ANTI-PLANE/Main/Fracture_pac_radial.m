clear all;clear all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KII = 0.5;
G   = 0;
l = 0.1;
eta = 0.9;
x1 = 0.1*l;
x2 = 0.4*l;
npt = 100;
%Radius = 1;
Radius = 20*l;

Ls = l/sqrt(2);
Lt = sqrt(eta+1)*l;
C  = G;
C1 = Ls*Ls;
C2 = Lt*Lt/2;
% Definizione della geometria
crv = nrbline([0 0],[Radius 0]);
srf = nrbrevolve(crv,[0 0 0],[0 0 Radius],pi);
srf = nrbdegelev(srf,[1 2]);
% 
% new_knots = linspace (0, 1, 10);  % Mettere dispari!!!!
% new_knots = new_knots (2:end-1);
% 
% srf_r = nrbkntins(srf, {new_knots, new_knots});
% problem_data.geo_name = srf_r;

problem_data.geo_name = srf;
% ----------------------------------------------------------------------- %
% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [1 2 4] ; %
problem_data.drchlt_sides = 4 ; % -> arco(far field solution)
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) C.*ones(size(x)) ;           % Laplace
problem_data.d_diff  = @(x, y) C1.*ones(size(x)) ;          % Bilaplace
problem_data.copp = C2 ;                                    % Bordi
% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt2;
 
% Derivate prime delle normali.
problem_data.dnorm_x_x  = @test_circular_plate_couple_dnx_x;
problem_data.dnorm_x_y  = @test_circular_plate_couple_dnx_y;
problem_data.dnorm_y_x  = @test_circular_plate_couple_dny_x;
problem_data.dnorm_y_y  = @test_circular_plate_couple_dny_y;
%
%problem_data.dnorm_x_xx = @test_circular_plate_couple_dnx_xx;
%problem_data.dnorm_x_xy = @test_circular_plate_couple_dnx_xy;
%problem_data.dnorm_x_yy = @test_circular_plate_couple_dnx_yy;
%
%problem_data.dnorm_y_xx = @test_circular_plate_couple_dny_xx;
%problem_data.dnorm_y_xy = @test_circular_plate_couple_dny_xy;
%problem_data.dnorm_y_yy = @test_circular_plate_couple_dny_yy;

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [10 10];     % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% method_data.degree     = [srf_r.dim-1 srf_r.dim-1];       % Degree of the splines
% method_data.regularity = [srf_r.dim-2 srf_r.dim-2];       % Regularity of the splines
% method_data.nsub       = [1 1];                           % Number of subdivisions
% method_data.nquad      = [4 4];                           % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_fracture_antiPlane_radial (problem_data, method_data);
%solve_fracture_antiPlane_NURBS_GRAD (problem_data, method_data);
% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, npt), linspace(0, 1, npt)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% Solution
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure, surf (X, Y, eu)
title ('Numerical solution','fontsize',20)
grid on
axis equal
axis on
axis tight

[eug, Fg] = sp_eval (u, space, geometry, vtk_pts, 'gradient');
[X, Y]  = deal (squeeze(Fg(1,:,:)), squeeze(Fg(2,:,:)));
gradx = squeeze(eug(1,:,:));
grady = squeeze(eug(2,:,:));

% figure, surf (X, Y, gradx)
% title ('Gradient: first component','fontsize',20)
% grid on
% axis equal
% axis on
% axis tight
% 
% figure, surf (X, Y, grady)
% title ('Gradient: second component','fontsize',20)
% grid on
% axis equal
% axis on
% axis tight

% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
n = size(X,2);
leap_SG_COMP(:,1) = [X(1,n:-1:1) X(npt,:)] ;
leap_SG_COMP(:,2) = [eu(1,n:-1:1) eu(npt,:)] ;

% Soluzione 1 Gradiente
anal_sol = (KII/G).*sqrt(-leap_SG_COMP(:,1).*(2/pi)) ;
%
% Soluzione 2 Gradiente
CIII = 2.8;
rr = linspace(0,1);
s = CIII.*rr.^(3/2).*( 3./ (16.*(Ls/Lt).^2-3) + 1 ) ;
%
% Lip
d = Lt*sqrt(3)/(16*sqrt(2)) *sqrt((Lt/Ls)^4 -32*(Lt/Ls)^2 +128) ;
disp(d)
% ----------------------------------------------------------------------- %
% Plot solution in any direction
j=1;
np = 100 ;
theta = zeros(1,np+1);
w1 = zeros(1,np+1);
w2 = zeros(1,np+1);
dt = 1/np;
for pn=0:dt:1
    theta(j) = 4*atand(pn);
    
    if pn==0.0
        a1 = pn ;
        a2 = pn ;
        b1 = x1/Radius;
        b2 = x2/Radius;
    elseif pn==1.0
        a1 = 1 ;
        a2 = 1 ;
        b1 = x1/Radius ;
        b2 = x2/Radius ;
    else
        a1 = pn/Radius*x1 ;
        a2 = pn/Radius*x2 ;
        b1 = 1/Radius*x1 ;
        b2 = 1/Radius*x2 ;
    end
    
    [uel1, p1] = sp_eval (u, space, geometry, {a1,b1});
    [uel2, p2] = sp_eval (u, space, geometry, {a2,b2});
    
    w1(j) = uel1 ;
    w2(j) = uel2 ;
    j= j+1;
end
figure, plot(theta,w1*G*sqrt(pi/(2*l))/KII)
title ('FIGURE 11 a) Normalized angular variation of the out-of-plane displacement field w at distance 0.1l','fontsize',15)
xlabel('\theta','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([0 180 0 0.4])
grid on

figure, plot(theta,w2*G*sqrt(pi/(2*l))/KII)
title ('FIGURE 11 b) Normalized angular variation of the out-of-plane displacement field w at distance 0.4l','fontsize',15)
xlabel('\theta','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([0 180 0 0.8])
grid on

% Valutazione flesso sulla lip
p=polyfit(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),9);
y=polyval(p,leap_SG_COMP(npt:end,1));
figure, plot(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),leap_SG_COMP(npt:end,1),y);

dp=polyder(p) ;
ddp=polyder(dp);
flex=roots(ddp);
disp(abs(max(flex)))
% ----------------------------------------------------------------------- %
figure
plot(leap_SG_COMP(:,1)./l,leap_SG_COMP(:,2).*(G*sqrt(pi/(2*l))/KII),'- r', 'linewidth',2 );
title ('FIGURE 9 a) Variation of crack face sliding displacement w along the crack face','fontsize',15)
xlabel('x_{1}/l','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([-5 0 0 2.5])
grid on

figure
plot(log(-leap_SG_COMP(:,1)./l), log(leap_SG_COMP(:,2).*(G*sqrt(pi/(2*l))/KII)),'- r', 'linewidth',2 );
title ('FIGURE 9 b) Variation of crack face sliding displacement w along the crack face in logarithmic scales','fontsize',15)
xlabel('log_{10}(x_{1}/l)','fontsize',20)
ylabel('log_{10}[wG(\pi/2l)^{1/2}/K_{III}]','fontsize',20)
axis([-3 1 -3 1])
grid on


% Plot Solution
%figure
%plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),leap_SG_COMP(:,1),anal_sol,-r,s)
%plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),'- r', 'linewidth',2 );
%gtext('Lip Solution' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');
%title ('Solution on the Lip','fontsize',20)
%xlabel('lip x \in [-1, 1]','fontsize',15)
%ylabel('Transverse Displacement','fontsize',15)
%axis equal
%grid on

%PLOT IN MATLAB.
% % % Plottagio geometria e control points
% figure
% nrbctrlplot(srf)
% axis on
% grid off
% %print('geometry.png')
% 
% % Plottaggio raffinamento
% figure
% nrbkntplot(srf)
% axis on
% grid off
% 
%print('geometry_cpoint.png')
%figure,
%nrbctrlplot(srf_r)
%figure
%nrbkntplot(srf_r)
% ----------------------------------------------------------------------- %
% figure
% plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),'- b', 'linewidth',1 );
% hold on
% plot(-rr,s, '- r', 'linewidth',1);
% hold on
% plot([-d -d],[0 0.8],'- g', 'linewidth',1 );
% hold off
% 
% title ('Solution on the Lip','fontsize',20)
% xlabel('lip x \in [-1, 1]','fontsize',15)
% ylabel('Transverse Displacement','fontsize',15)
% axis equal
% grid on
% 
% leap_SG_COMP0_100 = leap_SG_COMP;
% save leap_SG_COMP0_100
% ----------------------------------------------------------------------- %
% disp(msh.nel)
% format long
% disp( leap_SG_COMP(50,:)/Ls );
% disp( leap_SG_COMP(56,:)/Ls );

% ----------------------------------------------------------------------- %
% figure
% plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),'- b', 'linewidth',2 );
% hold on
% plot(-r,s, '- r', 'linewidth',2);
% hold on
% plot([-d -d],[0 0.8],'- k', 'linewidth',2 );
% hold off
% 
% gtext( 'Limit zone of influence' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
% gtext( 'of the second gradient (Analytic)' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
% gtext( 'Numerical solution' , 'fontsize', 20, 'rotation', 0, 'color' , 'b');
% gtext( 'Analytical solution' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');
% 
% title ('Lip Solution (Zoom)','fontsize',20)
% xlabel('lip x \in [-1, 1]','fontsize',15)
% ylabel('Transverse Displacement','fontsize',15)
% axis equal
% grid on
% save leap_SG_COMP