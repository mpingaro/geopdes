clear all ;
close all ;
clc ;
%--- Straight Beam problem -----------------------------------------------%
w = 1;
L = 100;
s = 100;

deg_u = 2;
deg_v = 2;
deg_w = 2;

nel_u = 8;
nel_v = 1;
nel_w = 1;

t = L/s;
% Geometry
line = nrbline([0 0 0],[L 0 0]);
srf  = nrbextrude(line, [0 w 0]);
vol  = nrbextrude(srf, [0 0 t]);
%figure, nrbplot(vol,[1 1 1]);

% Degree elevation
init_deg = vol.order-1;
if init_deg(1) == deg_u
    deg_elev(1) = 0;
else
    deg_elev(1) = deg_u-init_deg(1);
end

if init_deg(2) == deg_v
    deg_elev(2) = 0;
else
    deg_elev(2) = deg_v-init_deg(2);
end

if init_deg(1) == deg_u
    deg_elev(3) = 0;
else
    deg_elev(3) = deg_w-init_deg(3);
end
vol_r = nrbdegelev(vol, deg_elev);

% Knots insertion
knots_u = linspace(0,1,nel_u+1);
knots_u = knots_u(2:end-1);
knots_v = linspace(0,1,nel_v+1);
knots_v = knots_v(2:end-1);
knots_w = linspace(0,1,nel_w+1);
knots_w = knots_w(2:end-1);
vol_r = nrbkntins(vol_r, {knots_u knots_v knots_w});

% Extract value
format long
cp = reshape(vol_r.coefs(1:3,:,:,:),3,[]);
knot{1} = vol_r.knots{1}(deg_u+1:end-deg_u);
knot{2} = vol_r.knots{2}(deg_v+1:end-deg_v);
knot{3} = vol_r.knots{3}(deg_w+1:end-deg_w);

% Force
F = 250/((deg_u+1)*s^3);
tcp = size(cp,2);
tcp3 = 3*tcp;

fid = fopen('cantilever.txt','w');

fprintf(fid, 'Total number of control points 3D %12.0f\n', tcp3);
fprintf(fid, 'Total number of control points %12.0f\n', tcp);
fprintf(fid, 'Force applicate %12.12f\n', F);

fprintf(fid,'Knots direction u\n');
fprintf(fid,'%9.6f',knot{1});
fprintf(fid,'\n');

fprintf(fid,'Knots direction v\n');
fprintf(fid,'%9.6f',knot{2});
fprintf(fid,'\n');

fprintf(fid,'Knots direction w\n');
fprintf(fid,'%9.6f',knot{3});
fprintf(fid,'\n');

fprintf(fid,'Control points direction x\n');
fprintf(fid,'%11.6f',cp(1,:));
fprintf(fid,'\n');
fprintf(fid,'Control points direction y\n');
fprintf(fid,'%11.6f',cp(2,:));
fprintf(fid,'\n');
fprintf(fid,'Control points direction z\n');
fprintf(fid,'%11.6f',cp(3,:));

fclose(fid);

% %--- Circular Beam problem -----------------------------------------------%
% 
% t = 1 ; 
% w = 1 ;
% R = 10 ;
% 
% crv = nrbline([R, 0, 0],[R+t, 0, 0]);
% srf_sup = nrbrevolve(crv, [0 0 0], [0 -1 0], pi/2); 
% srf = nrbextrude(srf_sup, [0 w 0]) ;
% 
% figure, nrbplot(srf,[20 20 20]);
% %-------------------------------------------------------------------------%
% 
% %--- Pinched Cylinder problem --------------------------------------------%
% clear all ;
% 
% t =   3.00 ; 
% w = 600.00 ;
% R = 300.00 ;
% 
% crv = nrbline([R, 0, 0],[R+t, 0, 0]);
% srf_sup = nrbrevolve(crv, [0 0 0], [0 -1 0], pi/2); 
% srf = nrbextrude(srf_sup, [0 w/2 0]) ;
% 
% figure, nrbplot(srf,[30 3 30]);
% %-------------------------------------------------------------------------%
% 
% %--- Full hemispherical shell problem ------------------------------------%
% clear all ;
% 
% R = 10.00 ;
% t =  0.04 ;
% crv_line = nrbline([R-t/2 0 0],[R+t/2 0 0]) ;
% srf_arc = nrbrevolve(crv_line,[0 0 0],[0 -1 0],pi/2) ;
% srf = nrbrevolve(srf_arc,[0 0 0],[0 0 1],pi/2) ;
% 
% figure, nrbplot(srf,[20 20 20]) ;
% %-------------------------------------------------------------------------%
% 
% %--- Scordelis-Lo roof problem -------------------------------------------%
% clear all ;
% clc ;
% 
% R = 25.00 ;
% w = 50.00 ;
% t =  0.25 ;
% 
% crv = nrbline([0, 0, R],[0, 0, R+t]);
% srf_sup = nrbrevolve(crv, [0 0 0], [0 1 0], 2*pi/9); 
% srf = nrbextrude(srf_sup, [0 w/2 0]) ;
%  
% figure, nrbplot(srf,[25 2 25]);
% %-------------------------------------------------------------------------%