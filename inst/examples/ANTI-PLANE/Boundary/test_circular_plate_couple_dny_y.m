function dny_y = test_circular_plate_couple_dny_y (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dny_y = y.^2 ./ ( x.^3 .*(1 + y.^2./x.^2) ) + 1./( x.*(1 + y.^2./x.^2).^(1/2) );
      %dny_y = cos(theta) .* ( x ./ ( x.^2 + y.^2 ) );
      %dny_y = cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* 1 ./ x);
    case 2
      dny_y = y.^2 ./ ( x.^3 .*(1 + y.^2./x.^2) ) + 1./( x.*(1 + y.^2./x.^2).^(1/2) );  
      %dny_y = cos(theta) .* ( x ./ ( x.^2 + y.^2 ) );
      %dny_y = cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* 1 ./ x);
    case 3
      dny_y = y.^2 ./ ( x.^3 .*(1 + y.^2./x.^2) ) + 1./( x.*(1 + y.^2./x.^2).^(1/2) );  
      %dny_y = cos(theta) .* ( x ./ ( x.^2 + y.^2 ) );
      %dny_y = cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* 1 ./ x);
    case 4
      dny_y = y.^2 ./ ( x.^3 .*(1 + y.^2./x.^2) ) + 1./( x.*(1 + y.^2./x.^2).^(1/2) );  
      %dny_y = cos(theta) .* ( x ./ ( x.^2 + y.^2 ) );  
      %dny_y = cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* 1 ./ x);
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  if x<0
      dny_y = -dny_y;
  elseif x == 0
      dny_y = 0.*x.*y;
  end

end