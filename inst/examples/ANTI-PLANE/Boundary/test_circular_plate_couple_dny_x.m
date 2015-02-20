function dny_x = test_circular_plate_couple_dny_x (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dny_x =  y.^3 ./( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - y./( x.^2 .*( 1 + y.^2./x.^2).^(1/2) ); 
      %dny_x = - cos(theta) .* ( y ./ (x.^2 + y.^2) ); 
      %dny_x = - cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* y ./ x.^2 );
    case 2
      dny_x =  y.^3 ./( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - y./( x.^2 .*( 1 + y.^2./x.^2).^(1/2) );  
      %dny_x = - cos(theta) .* ( y ./ (x.^2 + y.^2) );  
      %dny_x = - cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* y ./ x.^2 );
    case 3
      dny_x =  y.^3 ./( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - y./( x.^2 .*( 1 + y.^2./x.^2).^(1/2) );  
      %dny_x = - cos(theta) .* ( y ./ (x.^2 + y.^2) );  
      %dny_x = - cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* y ./ x.^2 );
    case 4
      dny_x =  y.^3 ./( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - y./( x.^2 .*( 1 + y.^2./x.^2).^(1/2) );  
      %dny_x = - cos(theta) .* ( y ./ (x.^2 + y.^2) );  
      %dny_x = - cos(theta) .* ( 1 ./ (1 + (y./x).^2) .* y ./ x.^2 );
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<0
      dny_x = -dny_x;
  elseif x == 0
      dny_x = 0.*x.*y;
  end

end