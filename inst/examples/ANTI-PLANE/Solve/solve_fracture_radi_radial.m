function [geometry, msh, space, u] = ...
              solve_fracture_radi_radial (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', false);
  
% Construct space structure
space  = sp_nurbs_2d (geometry.nurbs, msh);

% Assemble the matrices
stiff_mat = op_fracture_radi_tp (space, space, msh, c_diff, d_diff, e_diff);
rhs       = op_f_v_tp (space, msh, f);

% Apply Neumann boundary conditions
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  sp_side  = sp_eval_boundary_side (space, msh_side);

  x = squeeze (msh_side.geo_map(1,:,:));
  y = squeeze (msh_side.geo_map(2,:,:));
  gval = reshape (g (x, y, iside), msh_side.nqn, msh_side.nel);

  rhs(sp_side.dofs) = rhs(sp_side.dofs) + op_f_v (sp_side, msh_side, gval);
end 

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs_u] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides_u);
u(drchlt_dofs_u) = u_drchlt;
%
[r_drchlt, drchlt_dofs_r] = sp_drchlt_l2_proj_gradudotn (space, msh, geometry, u, g, drchlt_sides_r);
u(drchlt_dofs_r) = r_drchlt;

% Apply Dirichlet boundary conditions (edge 2)
[u_drchlt2, drchlt_dofs2] = sp_drchlt_l2_proj (space, msh, h, 2);
ind = 1;
drchlt_dofs2 = drchlt_dofs2(ind);
u_drchlt2 = u_drchlt2(ind);
u(drchlt_dofs2) = u_drchlt2;

drchlt_dofs = union(drchlt_dofs_u,drchlt_dofs_r);
drchlt_dofs = union(drchlt_dofs, drchlt_dofs2);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
%
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs_u)*u_drchlt...
    -stiff_mat(int_dofs, drchlt_dofs_r)*r_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end