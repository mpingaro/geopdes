function h = boundary_pacman_mixed_bc_h_drchlt_new (x,y,ind)
  KIII = 1; % Lo tengo fisso!!
  
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
       h = 2.*KIII.*sqrt(r./(2.*pi)).*sin(theta./2);
    case 2  
       h = 2.*KIII.*sqrt(r./(2.*pi)).*sin(theta./2);
    case 3
      h = 0.*x + 0.*y;
    case 4
      h = 2.*KIII.*sqrt(r./(2.*pi)).*sin(theta./2);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end