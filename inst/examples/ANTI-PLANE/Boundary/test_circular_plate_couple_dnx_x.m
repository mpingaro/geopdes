function dnx_x = test_circular_plate_couple_dnx_x (x, y, ind)

  %[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dnx_x = y.^2./( x.^3 .*(1 + y.^2./x.^2).^(3/2) );
      %dnx_x = sin(theta) .* ( y./ ( x.^2 + y.^2) );
      %dnx_x = sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (y ./ x.^2) );
    case 2
      dnx_x = y.^2./( x.^3 .*(1 + y.^2./x.^2).^(3/2) );
      %dnx_x = sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (y ./ x.^2) );
      %dnx_x = sin(theta) .* ( y./ ( x.^2 + y.^2) );
    case 3
      dnx_x = y.^2./( x.^3 .*(1 + y.^2./x.^2).^(3/2) );
      %dnx_x = sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (y ./ x.^2) );
      %dnx_x = sin(theta) .* ( y./ ( x.^2 + y.^2) );
    case 4
      dnx_x = y.^2./( x.^3 .*(1 + y.^2./x.^2).^(3/2) );
      %dnx_x = sin(theta) .* ( 1 ./ (1+ (y ./ x).^2) .* (y ./ x.^2) );
      %dnx_x = sin(theta) .* ( y./ ( x.^2 + y.^2) );
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<=0
      dnx_x = -dnx_x;
  elseif x == 0
      dnx_x = 0.*x.*y;
  end
  
end