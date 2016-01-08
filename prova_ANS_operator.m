% By Marco Pingaro
clear;
clc;

% Initial Space
lin = nrbline([-1,-1],[1,-1]);
crv_0 = nrbextrude(lin, [0,2,0]);

%% Set of Gauss Points in the Parametric Domanin (0,1)--(0,1)
a = sqrt(3/5)/2;
b = sqrt(1/3)/2;

gp = {{0.5-a, 0.5-a};
      {0.5, 0.5-a};
      {0.5+a, 0.5-a};
      {0.5-a, 0.5};
      {0.5, 0.5};
      {0.5+a, 0.5};
      {0.5-a, 0.5+a};
      {0.5, 0.5+a};
      {0.5+a, 0.5+a}};  
  
%% First Local Space
crv_1 = nrbdegelev(crv_0, [0,1]);
tying_1 = {{0.5-b, 0.5-a};
           {0.5+b, 0.5-a};
           {0.5-b, 0.5};
           {0.5+b, 0.5};
           {0.5-b, 0.5+a};
           {0.5+b, 0.5+a}};
       
M_1 = zeros(6);
for i=1:6
    M_1(i,:) = nrbbasisfun(tying_1{i}, crv_1);
end
N_1 = zeros(9,6);
for i=1:9
    N_1(i,:) = nrbbasisfun(gp{i}, crv_1);
end
L_1 = N_1*inv(M_1);

%% Second Local Space
crv_2 = nrbdegelev(crv_0, [1,0]);
% Set of tying points
tying_2 = {{0.5-a, 0.5-b};
           {0.5, 0.5-b};
           {0.5+a, 0.5-b};
           {0.5-a, 0.5+b};
           {0.5, 0.5+b};
           {0.5+a, 0.5+b}};
M_2 = zeros(6);
for i=1:6
    M_2(i,:) = nrbbasisfun(tying_2{i}, crv_2);
end
N_2 = zeros(9,6);
for i=1:9
    N_2(i,:) = nrbbasisfun(gp{i}, crv_2);
end
L_2 = N_2*inv(M_2);

%% Third Local Space
crv_3 = nrbdegelev(crv_0, [0,0]);
% Set of tying points
tying_3 = {{0.5-b, 0.5-b};
           {0.5+b, 0.5-b};
           {0.5-b, 0.5+b};
           {0.5+b, 0.5+b}};
M_3 = zeros(4);
for i=1:4 
    M_3(i,:) = nrbbasisfun(tying_3{i}, crv_3);
end
N_3 = zeros(9,4);
for i=1:9
    N_3(i,:) = nrbbasisfun(gp{i}, crv_3);
end
L_3 = N_3*inv(M_3);

%% Computate da Josef
% for ANS: extrapolation matrices
% 1. for xixi and xizeta
M1 = [1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0 0 0;...
      1/2             1/2            0 0 0 0;...
      1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0 0 0;...
      0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20) 0 0;...
      0 0 1/2             1/2            0 0;...
      0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20) 0 0;...
      0 0 0 0 1/2+sqrt(9/20)  1/2-sqrt(9/20);...
      0 0 0 0 1/2             1/2           ;...
      0 0 0 0 1/2-sqrt(9/20)  1/2+sqrt(9/20)];
% 2. for etaeta and etazeta
M2 = [1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0 0;...
      0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20) 0;...
      0 0 1/2+sqrt(9/20) 0 0 1/2-sqrt(9/20);...
      1/2            0 0 1/2            0 0;...
      0 1/2            0 0 1/2            0;...
      0 0 1/2            0 0 1/2           ;...
      1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0 0;...
      0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20) 0;...
      0 0 1/2-sqrt(9/20) 0 0 1/2+sqrt(9/20)];
% 3. for xieta
M3 = 1/4*[14/5+sqrt(36/5) -4/5          -4/5          14/5-sqrt(36/5);...
          1+sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1-sqrt(9/5);...
          -4/5          14/5+sqrt(36/5) 14/5-sqrt(36/5) -4/5;...
          1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5);...
          1           1           1           1;...
          1-sqrt(9/5)   1+sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5);...
          -4/5          14/5-sqrt(36/5) 14/5+sqrt(36/5) -4/5;...
          1-sqrt(9/5)   1-sqrt(9/5)   1+sqrt(9/5)   1+sqrt(9/5);...
          14/5-sqrt(36/5) -4/5          -4/5          14/5+sqrt(36/5)];