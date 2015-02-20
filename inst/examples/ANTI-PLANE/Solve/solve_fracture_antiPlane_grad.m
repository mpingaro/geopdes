function [geometry, msh, space, u] = ...
              solve_fracture_antiPlane_grad (problem_data, method_data)

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
msh      = msh_2d (zeta, qn, qw, geometry,'der2', true);
  
% Construct space structure
space  = sp_nurbs_2d (geometry.nurbs, msh);

% Assemble the matrices
stiff_mat1 = op_fracture_antiPlaneGRAD_tp (space, space, msh, c_diff, d_diff);
rhs       = op_f_v_tp (space, msh, f);

% Apply Neumann boundary conditions

stiff_mat2 = zeros(space.ndof, space.ndof);
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  sp_side  = sp_eval_boundary_side (space, msh_side);
  

  x = squeeze (msh_side.geo_map(1,:,:));
  y = squeeze (msh_side.geo_map(2,:,:));
  gval = reshape (g (x, y, iside), msh_side.nqn, msh_side.nel);
  rhs(sp_side.dofs) = rhs(sp_side.dofs) + op_f_v (sp_side, msh_side, gval);
  
   if iside == 1
      rule_side = rule;
      rule_side{(1)} = [-1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', false);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_col = msh_evaluate_col(msh_bound,1);
      
      msh_col.normal = msh_side.normal;                  % Import of structure normal.
      msh_col.normal(1,:,:) = 0;
      
      sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true, 'hessian', true, 'der3', false);
      % Normal and Tangent derivatives 
      sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true,'hessian', true, 'der3', false);
     
      sp_boundary = sp_side_col;
      
      bordo = copp .* op_Bordo_Fracture(sp_boundary, sp_normal, msh_side);
      %bordo2 = copp .* op_Bordo_Fracture_grad(sp_boundary, sp_boundary, msh_side);
            
   elseif iside == 2
      rule_side = rule;
      rule_side{(1)} = [1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', true);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_col = msh_evaluate_col(msh_bound,msh.nel_dir(1));
     
      msh_col.normal = msh_side.normal;                  % Import of structure normal.
      msh_col.normal(1,:,:) = 0;
      
      sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true, 'hessian', true, 'der3', true);
      % Normal and Tangent derivatives  
      sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true,'hessian', true, 'der3', true);
      
      sp_boundary = sp_side_col;
      
      %bordo2 = copp .* op_Bordo_Fracture(sp_boundary, sp_normal, msh_side);
      bordo = copp .* op_Bordo_Fracture_grad(sp_boundary, sp_boundary, msh_side);
  
   elseif iside == 4
      % NORMAL DERIVATION CIRCULAR PART
      normal.deriv(1,1,:,:)  = reshape (dnorm_x_x (x, y, iside), msh_side.nqn, msh_side.nel);
      normal.deriv(1,2,:,:)  = reshape (dnorm_x_y (x, y, iside), msh_side.nqn, msh_side.nel);
      normal.deriv(2,1,:,:)  = reshape (dnorm_y_x (x, y, iside), msh_side.nqn, msh_side.nel);
      normal.deriv(2,2,:,:)  = reshape (dnorm_y_y (x, y, iside), msh_side.nqn, msh_side.nel);
    
      %normal.deriv2(1,1,:,:) = reshape (dnorm_x_xx (x, y, iside), msh_side.nqn, msh_side.nel);
      %normal.deriv2(1,2,:,:) = reshape (dnorm_x_xy (x, y, iside), msh_side.nqn, msh_side.nel);
      %normal.deriv2(1,3,:,:) = reshape (dnorm_x_yy (x, y, iside), msh_side.nqn, msh_side.nel);
    
      %normal.deriv2(2,1,:,:) = reshape (dnorm_y_xx (x, y, iside), msh_side.nqn, msh_side.nel);
      %normal.deriv2(2,2,:,:) = reshape (dnorm_y_xy (x, y, iside), msh_side.nqn, msh_side.nel);
      %normal.deriv2(2,3,:,:) = reshape (dnorm_y_yy (x, y, iside), msh_side.nqn, msh_side.nel); 
        
      rule_side = rule;
      rule_side{(2)} = [1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', false);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_row = msh_evaluate_row(msh_bound,msh.nel_dir(2));
      msh_row.normal = msh_side.normal;                  % Import of structure normal.
      
      sp_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true, 'hessian', true, 'der3', false);
      % Normal and Tangent derivatives
      sp_normal = sp_normal_tang_circ(sp_side_row, msh_row, normal, 'gradient', true,'hessian', true, 'der3', false);
      
      sp_boundary = sp_side_row;
      
      bordo = copp .* op_Bordo_Fracture(sp_boundary, sp_normal, msh_side);
   end
   
   %stiff_mat2 = stiff_mat2 + copp .* op_Bordo_Fracture(sp_boundary, sp_normal, msh_side);
   %stiff_mat2 = stiff_mat2 + copp .* op_Bordo_Fracture_grad(sp_boundary, sp_boundary, msh_side);
   stiff_mat2 = stiff_mat2 + bordo;
   
end

stiff_mat = stiff_mat1 - stiff_mat2;

%   
% u = zeros (space.ndof, 1);

% % Apply Dirichlet boundary conditions (edge 1)
% [u_drchlt1, drchlt_dofs1] = sp_drchlt_l2_proj (space, msh, h, 1);
% 
% iniz = 1;
% fine = size(drchlt_dofs1,1);
% 
% drchlt_dofs1 = drchlt_dofs1(iniz : fine);
% u_drchlt1 = u_drchlt1(iniz : fine);
% 
% % Apply Dirichlet boundary conditions (edge 2)
% [u_drchlt2, drchlt_dofs2] = sp_drchlt_l2_proj (space, msh, h, 2);
% 
% ind = 1;
% drchlt_dofs2 = drchlt_dofs2(ind);
% u_drchlt2 = u_drchlt2(ind);
% 
% % Apply Dirichlet boundary conditions (edges 4)
% [u_drchlt4, drchlt_dofs4] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
% 
% % Aggiungo bordi
% drchlt_dofs_u = [drchlt_dofs1;drchlt_dofs2;drchlt_dofs4];
% u_drchlt = [u_drchlt1; u_drchlt2; u_drchlt4];
% u(drchlt_dofs_u) = u_drchlt;

% [r_drchlt, drchlt_dofs_r] = sp_drchlt_l2_proj_gradudotn (space, msh, geometry, u_drchlt, hr, drchlt_sides_r);
% u(drchlt_dofs_r) = r_drchlt;
% %
% drchlt_dofs = union(drchlt_dofs_u,drchlt_dofs_r);
% int_dofs = setdiff (1:space.ndof, drchlt_dofs);
% %
% rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs_u)*u_drchlt...
%     -stiff_mat(int_dofs, drchlt_dofs_r)*r_drchlt;
% 
% % Solve the linear system
% u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);


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