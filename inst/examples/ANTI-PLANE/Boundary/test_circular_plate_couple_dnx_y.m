function dnx_y = test_circular_plate_couple_dnx_y (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dnx_y = - y./ ( x.^2 .*( 1 + y.^2./x.^2 ).^(3/2) ); 
      %dnx_y = -sin(theta) .* x ./(x.^2 + y.^2);   
      %dnx_y = -sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (1 ./ x) );
    case 2
      dnx_y = - y./ ( x.^2 .*( 1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_y = -sin(theta) .* x ./(x.^2 + y.^2);  
      %dnx_y = -sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (1 ./ x) );
    case 3
      dnx_y = - y./ ( x.^2 .*( 1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_y = -sin(theta) .* x ./(x.^2 + y.^2);  
      %dnx_y = -sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (1 ./ x) );
    case 4
      dnx_y = - y./ ( x.^2 .*( 1 + y.^2./x.^2 ).^(3/2) );  
      %dnx_y = -sin(theta) .* x ./(x.^2 + y.^2);  
      %dnx_y = -sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (1 ./ x) );
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<0
      dnx_y = -dnx_y;
  elseif x == 0
      dnx_y = 0.*x.*y;
  end

end