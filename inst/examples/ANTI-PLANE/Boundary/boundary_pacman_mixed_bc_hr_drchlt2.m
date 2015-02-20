function hr = boundary_pacman_mixed_bc_hr_drchlt2 (x,y,ind)
  KIII = 0.5; % Fissato
  G    = 1;   % Fissato 
  [theta, r] = cart2pol (x,y);
  switch (ind)
    case 1
      hr = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);% - 0.5.*(KIII/G).*sqrt(2./pi).*sin(theta./2)./sqrt(r);
    case 2
      hr = 0.*x+0.*y;
    case 3
      hr = 0.*x+0.*y;
    case 4
      hr = (KIII/G).*sqrt(2.*r./pi).*sin(theta./2);% - 0.5.*(KIII/G).*sqrt(2./pi).*sin(theta./2)./sqrt(r);
      otherwise
          error ('h_drchlt: unknown reference number')
  end

end