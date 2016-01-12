function h = boundary_plate_g_drchlt (x,y,ind)
  
  switch (ind)
    case 1
      h = 0*x.*y;
    case 2
      h = pi/4+0.*x.*y;
    case 3
      h = 0.*x.*y;
    case 4
      h = 0.*x.*y;
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end