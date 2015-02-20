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
nqn_dir = p + 3;

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

% Compute some quantities using the tensor product structure
nqn  = prod (nqn_dir);
ndof = prod (ndof_dir);
nsh  = prod (nsh_dir);
nel  = prod (nel_dir);

% Start the loop to compute the stiffness matrix
A = sparse (ndof, ndof);
rhs = zeros (ndof, 1);
for elem = 1:nel
  [iel, jel] = ind2sub (nel_dir, elem);
% Connectivity
  conn_u = num{1}(:,iel);
  conn_v = num{2}(:,jel);
  conn_u_rep = repmat (conn_u, 1, nsh_dir(2));
  conn_v_rep = repmat (conn_v', nsh_dir(1), 1);
  conn = sub2ind (ndof_dir, conn_u_rep, conn_v_rep);
  conn = conn(:); % Store it as a vector

% % The same thing could be done using meshgrid
%  [conn_v_rep, conn_u_rep] = meshgrid (num{1}(:,iel), num{2}(:,jel));
%   conn = sub2ind (ndof_dir, conn_u_rep, conn_v_rep);
%   conn = conn(:); % Store it as a vector


% Quadrature weights (using a loop)
  quad_weights = zeros (1, nqn);
  for ind_qn = 1:nqn
    [iqn, jqn] = ind2sub (nqn_dir, ind_qn);
    quad_weights(ind_qn) = qw{1}(iqn, iel) * qw{2}(jqn,jel);
  end

% Quadrature weights. The computation is simpler using kron
  quad_weights = kron (qw{2}(:,jel), qw{1}(:,iel))';


% Shape functions values (using two loops)
  shp_u = Nfun{1}(:,:,iel);
  shp_v = Nfun{2}(:,:,jel);
  shape_functions = zeros (nqn, nsh);
  for ind_sh = 1:nsh
    [ish, jsh] = ind2sub (nsh_dir, ind_sh);
    for ind_qn = 1:nqn
      [iqn,jqn] = ind2sub (nqn_dir, ind_qn);
      shape_functions(ind_qn,ind_sh) = shp_u(iqn,ish) * shp_v(jqn,jsh);
    end
  end

% Shape function values (using kron)
  shape_functions = kron (shp_v, shp_u);


% Shape functions derivatives
  shg_u = Nder{1}(:,:,iel);
  shg_v = Nder{2}(:,:,jel);

  shape_function_gradients = zeros (2, nqn, nsh);
  shape_function_gradients(1,:,:) = kron (shp_v, shg_u);
  shape_function_gradients(2,:,:) = kron (shg_v, shp_u);


% Evaluate the parameterization at the quadrature points
  geo_map = zeros (2, nqn);
  for ind_qn = 1:nqn
    for ind_sh = 1:nsh
      [ish, jsh] = ind2sub (nsh_dir, ind_sh);
      geo_map(:,ind_qn) = geo_map(:,ind_qn) + shape_functions(ind_qn,ind_sh) * ...
        B(:, conn_u(ish), conn_v(jsh));
    end
  end


% Derivatives of the parameterization with respect to parametric coordinates
  geo_map_jac = zeros (2, 2, nqn);
  for ind_qn = 1:nqn
    for ind_sh = 1:nsh
      [ish, jsh] = ind2sub (nsh_dir, ind_sh);
      for idim = 1:2
        for ider = 1:2
          geo_map_jac(idim,ider,ind_qn) = geo_map_jac(idim,ider,ind_qn) + ...
            shape_function_gradients(ider,ind_qn,ind_sh) * ...
            B(idim, conn_u(ish), conn_v(jsh));
        end
      end
    end
  end


% Determinant of the Jacobian matrix
  jacdet = reshape (geo_map_jac(1,1,:) .* geo_map_jac(2,2,:) - ...
           geo_map_jac(1,2,:) .* geo_map_jac(2,1,:), 1, nqn);


% Passing to the physical domain
  for ind_qn = 1:nqn
    JinvT = inv (geo_map_jac(:,:,ind_qn))';
    for ind_sh = 1:nsh
      shape_function_gradients(:,ind_qn,ind_sh) = ...
          JinvT * shape_function_gradients(:,ind_qn,ind_sh);
    end
  end


% Assemble the stiffness matrix
  A_loc = zeros (nsh, nsh);
  for ind_sh = 1:nsh
    for jnd_sh = 1:nsh
      A_loc(ind_sh, jnd_sh) = sum (quad_weights .* jacdet .* ...
        sum(shape_function_gradients(:,:,ind_sh) .* ...
            shape_function_gradients(:,:,jnd_sh), 1), 2);
    end
  end
  A(conn,conn) = A(conn,conn) + A_loc;

% Assemble the right-hand side
  rhs_loc = zeros (nsh, 1);
  f = ones (1, nqn);
  for ind_sh = 1:nsh
    rhs_loc(ind_sh) = sum (quad_weights .* jacdet .* ...
      shape_functions(:,ind_sh).' .* f, 2);
  end
  rhs(conn) = rhs(conn) + rhs_loc;
  
end

% Find the boundary dofs
drchlt_dofs = unique ...
  ([sub2ind(ndof_dir, 1:ndof_dir(1), ones(1,ndof_dir(1))), ...
   sub2ind(ndof_dir, 1:ndof_dir(1), ndof_dir(2)*ones(1,ndof_dir(1))), ...
   sub2ind(ndof_dir, ones(1,ndof_dir(2)), 1:ndof_dir(2)), ...
   sub2ind(ndof_dir, ndof_dir(1)*ones(1,ndof_dir(2)), 1:ndof_dir(2))]);


% Solve the linear system
u = zeros (ndof, 1);
int_dofs = setdiff (1:ndof, drchlt_dofs);
u(int_dofs) = A(int_dofs,int_dofs) \ rhs(int_dofs);



% Postprocessing.
% For each point, compute the parameterization (F), and the solution (eu)
nplot = [30 30];
F  = zeros (2, nplot(1), nplot(2));
eu = zeros (nplot(1), nplot(2));

% Start computing the univariate fields
for idim = 1:2
  x{idim} = linspace (0, 1, nplot(idim));
  iv{idim} = findspan (ndof_dir(idim)-1, p(idim), x{idim}, knt{idim});
  N{idim} = basisfun (iv{idim}, x{idim}, p(idim), knt{idim});
  num{idim} = numbasisfun (iv{idim}(1,:), x{idim}, p(idim), knt{idim}) + 1;
end

% Then, compute the parameterization and the solution for every point
for ipt = 1:nplot(1)
  for jpt = 1:nplot(2)
    conn_u = num{1}(ipt,:);
    conn_v = num{2}(jpt,:);
    [conn_v_rep2, conn_u_rep2] = meshgrid (conn_v, conn_u);
    conn_u_rep = repmat (conn_u', 1, nsh_dir(2));
    conn_v_rep = repmat (conn_v, nsh_dir(1), 1);
    conn_plot = sub2ind (ndof_dir, conn_u_rep(:), conn_v_rep(:));

    shp_u = N{1}(ipt,:);
    shp_v = N{2}(jpt,:);
    shape_functions = kron (shp_v, shp_u);

% Evaluate the parameterization at the point
    for ind_sh = 1:nsh
      [ish, jsh] = ind2sub (nsh_dir, ind_sh);
      F(:,ipt,jpt) = F(:,ipt,jpt) + shape_functions(ind_sh) * ...
        B(:, conn_u(ish), conn_v(jsh));
    end
    
% Evaluate the computed solution at the point
    for ind_sh = 1:nsh
      eu(ipt,jpt) = eu(ipt,jpt) + shape_functions(ind_sh) * ...
        u(conn_plot(ind_sh));
    end
    
  end
end

surf(squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)
title ('Computed solution')
