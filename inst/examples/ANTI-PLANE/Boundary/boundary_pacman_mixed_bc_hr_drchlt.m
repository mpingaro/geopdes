function hr = boundary_pacman_mixed_bc_hr_drchlt (x,y,ind)
  KIII = 0.5;  % Fissato
  mu    = 1;   % Fissato 
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      hr = 0.*x+0.*y;
    case 2
      hr = 0.*x+0.*y;
    case 3
      hr = 0.*x+0.*y;
    case 4
      hr = (KIII/mu).*sqrt(2./pi).*sin(theta./2)./sqrt(r)./2;
      otherwise
          error ('h_drchlt: unknown reference number')
  end

end