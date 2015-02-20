function dny_yy = test_circular_plate_couple_dny_yy (x, y, ind)
%[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dny_yy = 3.*y.^3 ./ ( x.^5 .*(1 + y.^2./x.^2).^(5/2) ) - 3.*y ./ ( x.^3 .*(1 + y.^2./x.^2).^(3/2) ); 
      %dny_yy = -sin(theta) .* (x.^2 ./ (x.^2 + y.^2).^2) - cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );
      %dny_yy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 2
      dny_yy = 3.*y.^3 ./ ( x.^5 .*(1 + y.^2./x.^2).^(5/2) ) - 3.*y ./ ( x.^3 .*(1 + y.^2./x.^2).^(3/2) );   
      %dny_yy = -sin(theta) .* (x.^2 ./ (x.^2 + y.^2).^2) - cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );
      %dny_yy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 3
      dny_yy = 3.*y.^3 ./ ( x.^5 .*(1 + y.^2./x.^2).^(5/2) ) - 3.*y ./ ( x.^3 .*(1 + y.^2./x.^2).^(3/2) );   
      %dny_yy = -sin(theta) .* (x.^2 ./ (x.^2 + y.^2).^2) - cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );  
      %dny_yy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    case 4
      dny_yy = 3.*y.^3 ./ ( x.^5 .*(1 + y.^2./x.^2).^(5/2) ) - 3.*y ./ ( x.^3 .*(1 + y.^2./x.^2).^(3/2) );   
      %dny_yy = -sin(theta) .* (x.^2 ./ (x.^2 + y.^2).^2) - cos(theta) .* ( 2.*x.*y ./ (x.^2 + y.^2).^2 );  
      %dny_yy = -(cos(theta).*(x.^2 - y.^2))./(x.^2 + y.^2).^2;
    otherwise
      error ('g_nmnn: unknown reference number')
  end
  
  if x<0
      dny_yy = -dny_yy;
  elseif x == 0
      dny_yy = 0.*x.*y;
  end

end