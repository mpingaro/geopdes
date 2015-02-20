function dnx_yy = test_circular_plate_couple_dnx_yy (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dnx_yy = 3.*y.^2 ./( x.^4 .*(1 + y.^2./x.^2).^(5/2) ) - 1 ./ ( x.^2 .*( 1 + y.^2./x.^2).^(3/2) );  
      %dnx_yy = - cos(theta) .* ( x.^2 ./ (x.^2 + y.^2).^2 ) + sin(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2  );
      %dnx_yy = (2.*x.*y.*sin(theta))./(x.^2 + y.^2).^2;
    case 2
      dnx_yy = 3.*y.^2 ./( x.^4 .*(1 + y.^2./x.^2).^(5/2) ) - 1 ./ ( x.^2 .*( 1 + y.^2./x.^2).^(3/2) );  
      %dnx_yy = - cos(theta) .* ( x.^2 ./ (x.^2 + y.^2).^2 ) + sin(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2  );  
      %dnx_yy = (2.*x.*y.*sin(theta))./(x.^2 + y.^2).^2;
    case 3
      dnx_yy = 3.*y.^2 ./( x.^4 .*(1 + y.^2./x.^2).^(5/2) ) - 1 ./ ( x.^2 .*( 1 + y.^2./x.^2).^(3/2) );  
      %dnx_yy = - cos(theta) .* ( x.^2 ./ (x.^2 + y.^2).^2 ) + sin(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2  );  
      %dnx_yy = (2.*x.*y.*sin(theta))./(x.^2 + y.^2).^2;
    case 4
      dnx_yy = 3.*y.^2 ./( x.^4 .*(1 + y.^2./x.^2).^(5/2) ) - 1 ./ ( x.^2 .*( 1 + y.^2./x.^2).^(3/2) );  
      %dnx_yy = - cos(theta) .* ( x.^2 ./ (x.^2 + y.^2).^2 ) + sin(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2  );  
      %dnx_yy = (2.*x.*y.*sin(theta))./(x.^2 + y.^2).^2;
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<0
      dnx_yy = -dnx_yy;
  elseif x == 0
      dnx_yy = 0.*x.*y;
  end

end