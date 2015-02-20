% Compute the solution near the apex of the fracture 
% Cite Radi (81)
clear all; clc;

%% Costant
KIII = 0.5;
G    = 1.0;
l    = 0.1;
eta  = 0.9;
Radius = 20*l;

%% ------------------------------------------------------------------------
b = 3*( (1+eta)^2 - (1-eta/3)*sqrt(1+14/3*eta+9*eta^2) )^(1/3);
a =  abs ( sqrt( 1./(3-eta).*( 2 + 2^(7/3).*eta./b + b./(2^(1/3)*3) ) ) );

x = 0:5*l/10:5*l;
w = zeros(1,length(x));
for i = 1:length(x)
x1 = -x(i);
%H  = @(t) heaviside(t-1);
xmin = 0;
xmax = 1;
%Fr = @(t) atan( (t.^3.*sqrt(1-t.^2))./(1./(1+eta) - t.^2).^2 )./(t+t);
%R  = (1/pi)*integral(Fr,xmin,xmax );

%Fu = @(t) (exp( sqrt(2).*t.*x1./l ) -1 ).*( ( 1/(1+eta) - t.^2 ).^2 + t.^3.*sqrt(t.^2-1).*heaviside(t-1) ).*exp(-R);
%Fd = @(t) t.*sqrt(t).*( t.*sqrt(1+eta) -a ).*( (3-eta).*t.^4 - (4.*a^2-1)/(a.^4.*(1+eta)).*t.^2 + 1/(a.^2.*(1+eta).^2) );
%F = @(t) ((exp( sqrt(2).*t.*x1./l ) -1 ).*( ( 1/(1+eta) - t.^2 ).^2 + t.^3.*sqrt(t.^2-1).*heaviside(t-1) ).*exp(-R))./...
%    (t.*sqrt(t).*( t.*sqrt(1+eta) -a ).*( (3-eta).*t.^4 - (4.*a^2-1)/(a.^4.*(1+eta)).*t.^2 + 1/(a.^2.*(1+eta).^2) ));

F = @(t) ((exp( sqrt(2).*t.*x1./l ) - 1 ).*( ( 1/(1+eta) - t.^2 ).^2 + t.^3.*sqrt(t.^2-1).*heaviside(t-1) ).*exp( -(1./pi).*integral( @(s) atan( (s.^3.*sqrt(1-s.^2))./(1./(1+eta) - s.^2).^2 )./(s+t),xmin,xmax,'ArrayValued',true ) ) )./...
    (t.*sqrt(t).*( t.*sqrt(1+eta) -a ).*( (3-eta).*t.^4 - (4.*a^2-1)/(a.^4.*(1+eta)).*t.^2 + 1/(a.^2.*(1+eta).^2) ));


tmin = 0;
tmax = Inf;
w(i) = KIII*sqrt( l*(3-eta) )./( 2^(5/4)*pi*G ).*integral(F,tmin,tmax,'ArrayValued',true);
end
%% ------------------------------------------------------------------------
figure, plot(-x/l,w*G*sqrt( pi/(2*l) )/KIII,'b-.','LineWidth',1.4);
hold on
wf = KIII/G*sqrt(2*x/pi);
plot(-x/l, wf*G*sqrt(pi/(2*l))/KIII,'r:','LineWidth',1.4);
hold on
ws = ( KIII*x.^(3/2)/( G*l*sqrt(2*pi*(1+eta)*(3-eta)) ) ) * ( (1+eta)*(5/3+eta) );
plot(-x/l, ws*G*sqrt(pi/(2*l))/KIII,'k','LineWidth',1);
hold off
xlim([-5 0])
ylim([0 2.5])