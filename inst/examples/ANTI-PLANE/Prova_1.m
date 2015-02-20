clear all; close all; clc;
p = 3; 
knt = [0 0 0 1/5 2/5 3/5 4/5 1 1 1 1]; 
ncp = numel(knt)-(p+1); 
u = linspace (0, 1, 1000);  

iv = findspan (ncp-1, p, u, knt);  

N = basisfun (iv, u, p, knt);

num = numbasisfun(iv,u,p,knt);
num = num +1;

%Nplot = zeros(numel(u), ncp);
%for ii = 1:numel(u)
%    Nplot(ii,num(ii,:)) = N(ii,:);
%end
%plot(u, Nplot)


B = zeros(2,ncp);
B(1,:) = 1:ncp; B(2,:) = rand(1,ncp);
C = zeros(2, numel(u));

for ii = 1:numel(u)+1
    for ifun = 1: p+1
        C(:,ii) = C(:,ii) + N(ii,ifun) * B(:,num(ii,ifun));
    end
end

plot(C(1,:),C(2,:),'b-', B(1,:), B(2,:),'k--o')