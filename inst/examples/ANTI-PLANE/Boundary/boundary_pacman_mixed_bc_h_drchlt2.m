function h = boundary_pacman_mixed_bc_h_drchlt2 (x,y,ind)
  KIII = 0.5; % Fissato
  G    = 1;   % Fissato 
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      h = 0.*x+0.*y;
    case 2
      h = 0.*x+0.*y;
    case 3
      h = 0.*x+0.*y;
    case 4
      h = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);
    otherwise
      error ('h_drchlt: unknown reference number')
  end

end