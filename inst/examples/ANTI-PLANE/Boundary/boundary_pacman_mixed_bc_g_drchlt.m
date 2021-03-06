function h = boundary_pacman_mixed_bc_g_drchlt (x,y,ind)
  KIII = 0.5; 
  mu    = 1;
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      h = ( (KIII/mu).*sqrt(2./pi).*sin(theta./2)./sqrt(r)./2 );
    case 2
      h = 0.*x.*y;
    case 3
      h = 0.*x.*y;
    case 4
      h = ( (KIII/mu).*sqrt(2./pi).*sin(theta./2)./sqrt(r)./2 );
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end


