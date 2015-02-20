function h = boundary_pacman_mixed_bc_h_drchlt (x,y,ind)
  KIII = 0.5; 
  G    = 1;
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      h = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);
    case 2
      h = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);
    case 3
      h = 0.*x.*y;
    case 4
      h = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end