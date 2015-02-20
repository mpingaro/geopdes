function dny_xy = test_circular_plate_couple_dny_xy (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dny_xy =  3.*y.^4 ./ ( x.^6 .*(1 + y.^2./x.^2).^(5/2) ) + 4.*y.^2 ./ ( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - 1 ./ ( x.^2 .*(1 + y.^2./x.^2).^(1/2) );
      %dny_xy = sin(theta) .* ( x.*y ./ (x.^2 + y.^2).^2 ) - cos(theta) .* ( (x.^2 - y.^2)./(x.^2 + y.^2).^2 );  
      %dny_xy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 2
      dny_xy =  3.*y.^4 ./ ( x.^6 .*(1 + y.^2./x.^2).^(5/2) ) + 4.*y.^2 ./ ( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - 1 ./ ( x.^2 .*(1 + y.^2./x.^2).^(1/2) );  
      %dny_xy = sin(theta) .* ( x.*y ./ (x.^2 + y.^2).^2 ) - cos(theta) .* ( (x.^2 - y.^2)./(x.^2 + y.^2).^2 );
      %dny_xy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 3
      dny_xy =  3.*y.^4 ./ ( x.^6 .*(1 + y.^2./x.^2).^(5/2) ) + 4.*y.^2 ./ ( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - 1 ./ ( x.^2 .*(1 + y.^2./x.^2).^(1/2) );  
      %dny_xy = sin(theta) .* ( x.*y ./ (x.^2 + y.^2).^2 ) - cos(theta) .* ( (x.^2 - y.^2)./(x.^2 + y.^2).^2 );  
      %dny_xy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 4
      dny_xy =  3.*y.^4 ./ ( x.^6 .*(1 + y.^2./x.^2).^(5/2) ) + 4.*y.^2 ./ ( x.^4 .*(1 + y.^2./x.^2).^(3/2) ) - 1 ./ ( x.^2 .*(1 + y.^2./x.^2).^(1/2) );  
      %dny_xy = sin(theta) .* ( x.*y ./ (x.^2 + y.^2).^2 ) - cos(theta) .* ( (x.^2 - y.^2)./(x.^2 + y.^2).^2 );  
      %dny_xy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    otherwise
      error ('g_nmnn: unknown reference number')
  end

  if x<0
      dny_xy = -dny_xy;
  elseif x == 0
      dny_xy = 0.*x.*y;
  end
  
end