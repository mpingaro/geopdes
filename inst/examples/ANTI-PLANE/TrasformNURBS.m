den_sum    = sum(N);  % N = N*w ho già  moltiplicato per i pesi.

der_sumx   = sum(dN(1,:));
der_sumy   = sum(dN(2,:));
der2_sumx  = sum(ddN(1,:));
der2_sumxy = sum(ddN(2,:));
der2_sumy  = sum(ddN(3,:));
                
ddN(1,:) = ddN(1,:)/den_sum - (2*dN(1,:)*der_sumx ...
    + N*der2_sumx)/den_sum^2 + 2*N*der_sumx^2/den_sum^3;
ddN(2,:) = ddN(2,:)/den_sum - (dN(1,:)*der_sumy ...
    + dN(2,:)*der_sumx + N*der2_sumxy)/den_sum^2 ...
    + 2*N*der_sumx*der_sumy/den_sum^3;
ddN(3,:) = ddN(3,:)/den_sum - (2*dN(2,:)*der_sumy ...
    + N*der2_sumy)/den_sum^2 + 2*N*der_sumy^2/den_sum^3;
dN(1,:)  = dN(1,:)/den_sum - N*der_sumx/den_sum^2;
dN(2,:)  = dN(2,:)/den_sum - N*der_sumy/den_sum^2;

N = N/den_sum;
                
X = reshape(Xl,nnodl,1);
Y = reshape(Yl,nnodl,1);

% compute derivatives with respect to physical coordinates
dxdxi = [dN(1,:)*X dN(2,:)*X;
         dN(1,:)*Y dN(2,:)*Y];

d2xdxi2 = [ddN(1,:)*X ddN(2,:)*X ddN(3,:)*X;
           ddN(1,:)*Y ddN(2,:)*Y ddN(3,:)*Y];
                
dxdxi2 = [dxdxi(1,1)^2            dxdxi(1,1)*dxdxi(1,2)                       dxdxi(1,2)^2;
 2*dxdxi(1,1)*dxdxi(2,1) dxdxi(1,1)*dxdxi(2,2)+dxdxi(1,2)*dxdxi(2,1) 2*dxdxi(1,2)*dxdxi(2,2);
     dxdxi(2,1)^2            dxdxi(2,1)*dxdxi(2,2)                       dxdxi(2,2)^2];
                
dN  = dxdxi'\dN;
ddN = dxdxi2'\(ddN - d2xdxi2'*dN);