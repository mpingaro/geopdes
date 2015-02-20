function dny_xx = test_circular_plate_couple_dny_xx (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dny_xx = 3.*y.^5 ./ ( x.^7 .*( 1 + y.^2./x.^2).^(5/2) ) - 5.*y.^3./ ( x.^5 .*( 1 + y.^2./x.^2).^(3/2) ) + 2.*y./( x.^3.*(1+y.^2./x.^2).^(1/2) ); 
      %dny_xx = -sin(theta) .* ( y.^2 ./ (x.^2 + y.^2).^2 ) + cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 ); 
      %dny_xx = (2.*x.*y.*cos(theta))./(x.^2 + y.^2).^2;
    case 2
      dny_xx = 3.*y.^5 ./ ( x.^7 .*( 1 + y.^2./x.^2).^(5/2) ) - 5.*y.^3./ ( x.^5 .*( 1 + y.^2./x.^2).^(3/2) ) + 2.*y./( x.^3.*(1+y.^2./x.^2).^(1/2) );  
      %dny_xx = -sin(theta) .* ( y.^2 ./ (x.^2 + y.^2).^2 ) + cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );  
      %dny_xx = (2.*x.*y.*cos(theta))./(x.^2 + y.^2).^2;
    case 3
      dny_xx = 3.*y.^5 ./ ( x.^7 .*( 1 + y.^2./x.^2).^(5/2) ) - 5.*y.^3./ ( x.^5 .*( 1 + y.^2./x.^2).^(3/2) ) + 2.*y./( x.^3.*(1+y.^2./x.^2).^(1/2) );  
      %dny_xx = -sin(theta) .* ( y.^2 ./ (x.^2 + y.^2).^2 ) + cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );  
      %dny_xx = (2.*x.*y.*cos(theta))./(x.^2 + y.^2).^2;
    case 4
      dny_xx = 3.*y.^5 ./ ( x.^7 .*( 1 + y.^2./x.^2).^(5/2) ) - 5.*y.^3./ ( x.^5 .*( 1 + y.^2./x.^2).^(3/2) ) + 2.*y./( x.^3.*(1+y.^2./x.^2).^(1/2) );  
      %dny_xx = -sin(theta) .* ( y.^2 ./ (x.^2 + y.^2).^2 ) + cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );  
      %dny_xx = (2.*x.*y.*cos(theta))./(x.^2 + y.^2).^2;
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<0
      dny_xx = -dny_xx;
  elseif x == 0
      dny_xx = 0.*x.*y;
  end

end