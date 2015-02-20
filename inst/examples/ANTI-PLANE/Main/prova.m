clear all; clc;

%syms x t;
eta = 0.9;

%f = simplify( atan( (t^3*sqrt(1-t^2))/(1/(1+eta) - t^2)^2 )/(x+t) );
%Rx = int(atan( (t^3*sqrt(1-t^2))/(1/(1+eta) - t^2)^2 )/(x+t),t,0,1,'PrincipalValue', true);
%f=int(exp(-Rx),x, 0,Inf);
%Rx = @(x) integral(@(t) atan( (t.^3.*sqrt(1-t.^2))./(1./(1+eta) - t.^2).^2 )./(x+t),0,1);

%Rx = @(x) integral(@(t) x + t,0,1,'ArrayValued',true);
%R  = integral(@(x) x*integral( @(t) x + t, 0,1 ),0,2,'ArrayValued', true);

R = integral( @(x) exp( -(1./pi).*integral( @(t) atan( (t.^3.*sqrt(1-t.^2))./(1./(1+eta) - t.^2).^2 )./(x+t), 0,1)),0,Inf,'ArrayValued', true);

% dx = 0.1;
% x = 0:dx:10^3;
% dt = 0.01;
% t = 0:dt:1;
% eta = 0.9;
% 
% R = zeros(size(x));
% Lx = length(R);
% Lt = length(t);
% 
% for ii=1:Lx
%     X = x(ii);
%     for jj=1:Lt
%         R(ii)=R(ii)+1/pi*dt*(atan(t(jj)^3*sqrt(1-t(jj)^2)/(1/(1+eta)-t(jj)^2)^2)/(t(jj)+X));
%     end
% end
% 
% b = 3*((1+eta)^2-(1-eta/3)*sqrt(1+14/3*eta+9*eta^2))^(1/3);
% a = sqrt(1/(3-eta)*(2+(2^(7/3)*eta)/b+b/(2^(1/3)*3)));
% elle = 0.1;
% 
% x1 = -5*elle:5*elle/100:0;
% tcrit = a/sqrt(1+eta); epsi = 0.001;
% dt1 = tcrit/100;
% dt2 = 1;
% tfin = 100;
% t1 = 0:dt1:tcrit-epsi;      L1 = length(t1);
% t2 = tcrit+epsi:dt2:tfin;   L2 = length(t2);
% 
% w = 0;
% for i=1:L1
%     w 
