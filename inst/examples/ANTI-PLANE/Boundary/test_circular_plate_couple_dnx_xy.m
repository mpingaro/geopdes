function dnx_xy = test_circular_plate_couple_dnx_xy (x, y, ind)

  %[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dnx_xy = - 3.*y.^3 ./ ( x.^5 .*( 1 + y.^2./x.^2).^(5/2) ) + 2.*y ./ ( x.^3 .*(1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_xy = cos(theta) .* ( x.*y ./(x.^2 + y.^2 ).^2 ) + sin(theta) .* ( (x.^2 - y.^2) ./ (x.^2 + y.^2).^2 ); 
      %dnx_xy = (sin(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 2
      dnx_xy = - 3.*y.^3 ./ ( x.^5 .*( 1 + y.^2./x.^2).^(5/2) ) + 2.*y ./ ( x.^3 .*(1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_xy = cos(theta) .* ( x.*y ./(x.^2 + y.^2 ).^2 ) + sin(theta) .* ( (x.^2 - y.^2) ./ (x.^2 + y.^2).^2 );  
      %dnx_xy = (sin(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 3
      dnx_xy = - 3.*y.^3 ./ ( x.^5 .*( 1 + y.^2./x.^2).^(5/2) ) + 2.*y ./ ( x.^3 .*(1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_xy = cos(theta) .* ( x.*y ./(x.^2 + y.^2 ).^2 ) + sin(theta) .* ( (x.^2 - y.^2) ./ (x.^2 + y.^2).^2 );  
      %dnx_xy = (sin(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 4
      dnx_xy = - 3.*y.^3 ./ ( x.^5 .*( 1 + y.^2./x.^2).^(5/2) ) + 2.*y ./ ( x.^3 .*(1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_xy = cos(theta) .* ( x.*y ./(x.^2 + y.^2 ).^2 ) + sin(theta) .* ( (x.^2 - y.^2) ./ (x.^2 + y.^2).^2 );  
      %dnx_xy = (sin(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    otherwise
      error ('g_nmnn: unknown reference number')
  end

  if x<0
      dnx_xy = -dnx_xy;
  elseif x == 0
      dnx_xy = 0.*x.*y;
  end
  
end