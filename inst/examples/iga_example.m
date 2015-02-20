p(1) = 2; p(2) = 3; % Degree
knt{1} = [0, 0, 0, 1/2, 1, 1, 1]; % Open knot vector
knt{2} = [0, 0, 0, 0, 1/3, 2/3, 1, 1, 1, 1]; % Open knot vector
% Control points to define the geometry
B(1,:,:) = [1, 1, 0.9, 0.5, 0.15, 0;
            1.25, 1.25, 1.15, 0.6, 0.2, 0;
            1.75 1.75, 1.6, 0.9, 0.3, 0;
            2, 2, 1.8, 1, 0.3, 0];
B(2,:,:) = fliplr (squeeze(B(1,:,:))/1.5);
        
% Number of control points and quadrature nodes
ndof_dir(1) = numel (knt{1}) - (p(1)+1);
ndof_dir(2) = numel (knt{2}) - (p(2)+1);
nqn_dir = p + 1;

% These are auxiliary functions to compute the quadrature points
rule = msh_gauss_nodes (nqn_dir);
[qn, qw] = msh_set_quad_nodes (knt, rule);

% Precompute everything in 1D
for idim = 1:2
  nel_dir(idim) = numel (unique (knt{idim})) - 1;
  iv{idim} = findspan (ndof_dir(idim)-1, p(idim), qn{idim}, knt{idim});
  N{idim} = basisfunder (iv{idim}, p(idim), qn{idim}, knt{idim}, 1);
  
% For the numbering it is enough to consider one quadrature point.
  num{idim} = numbasisfun (iv{idim}(1,:), qn{idim}(1,:), p(idim), knt{idim}) + 1;
  
% It is necessary to reshape the variable N
  nsh_dir(idim) = size (N{idim}, 3);
  Nfun{idim} = reshape (N{idim}(:,1,:), nqn_dir(idim), nel_dir(idim), nsh_dir(idim));
  Nder{idim} = reshape (N{idim}(:,2,:), nqn_dir(idim), nel_dir(idim), nsh_dir(idim));

% For convenience, I put the element as the last index  
  Nfun{idim} = permute (Nfun{idim}, [1, 3, 2]);
  Nder{idim} = permute (Nder{idim}, [1, 3, 2]);
  num{idim} = permute (num{idim}, [2, 1]);
end
