function h = boundary_pacman_mixed_bc_h_drchlt_radial (x,y,ind)
  KIII = 0.5; 
  mu    = 1;
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      h = 0.*x.*y;
    case 2
      h = 0.*x.*y;
    case 3
      h = 0.*x.*y;
    case 4
      h = (KIII/mu).*sqrt(2./pi).*sqrt(r).*sin(theta./2);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end