function [dnx_x, dnx_y, dny_x, dny_y] = test_circular_plate_couple (x, y, ind)
  %[theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      dnx_x = y^2/(x^3*(y^2/x^2 + 1)^(3/2));
      dnx_y = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      
      dny_x = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      dny_y = 1/(x*(y^2/x^2 + 1)^(3/2));
    case 2
      dnx_x = y^2/(x^3*(y^2/x^2 + 1)^(3/2));
      dnx_y = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      
      dny_x = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      dny_y = 1/(x*(y^2/x^2 + 1)^(3/2));
    case 3
      dnx_x = y^2/(x^3*(y^2/x^2 + 1)^(3/2));
      dnx_y = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      
      dny_x = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      dny_y = 1/(x*(y^2/x^2 + 1)^(3/2));
    case 4
      dnx_x = y^2/(x^3*(y^2/x^2 + 1)^(3/2));
      dnx_y = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      
      dny_x = -y/(x^2*(y^2/x^2 + 1)^(3/2));
      dny_y = 1/(x*(y^2/x^2 + 1)^(3/2));
    otherwise
      error ('g_nmnn: unknown reference number')
  end

end