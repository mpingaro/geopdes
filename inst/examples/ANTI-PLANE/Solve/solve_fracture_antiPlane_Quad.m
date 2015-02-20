function [geometry, msh, space, u] = ...
              solve_fracture_antiPlane_Quad (problem_data, method_data)

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
stiff_mat1 = op_fracture_antiPlane_tp (space, space, msh, c_diff, d_diff);
%rhs       = op_f_v_tp (space, msh, f);  % BULK LOAD !!!
rhs = zeros(space.ndof,1);

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
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', true);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_col = msh_evaluate_col(msh_bound,1);
      msh_col.normal = msh_side.normal;                  % Import of structure normal.
      
      sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true, 'hessian', true, 'der3', true);
      % Normal and Tangent derivatives 
      sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true,'hessian', true, 'der3', true);
      
      sp_boundary = sp_side_col;
      
   elseif iside == 2
      rule_side = rule;
      rule_side{(1)} = [1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', true);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_col = msh_evaluate_col(msh_bound,msh.nel_dir(1));
      msh_col.normal = msh_side.normal;                  % Import of structure normal.
      
      sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true, 'hessian', true, 'der3', true);
      % Normal and Tangent derivatives  
      sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true,'hessian', true, 'der3', true);
      
      sp_boundary = sp_side_col;
      
   elseif iside == 3
      rule_side = rule;
      rule_side{(2)} = [-1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', true);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_row = msh_evaluate_row(msh_bound,1);
      msh_row.normal = msh_side.normal;                  % Import of structure normal.
      
      sp_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true, 'hessian', true, 'der3', true);
      % Normal and Tangent derivatives
      sp_normal = sp_normal_tang(sp_side_row, msh_row,'gradient', true,'hessian', true, 'der3', true);
      
      sp_boundary = sp_side_row;
   
   elseif iside == 4
      rule_side = rule;
      rule_side{(2)} = [1; 2];
      [qn, qw] = msh_set_quad_nodes (zeta, rule_side);
      msh_bound = msh_2d (zeta, qn, qw, geometry,'der2', true, 'der3', true);
      space_side = sp_nurbs_2d (geometry.nurbs, msh_bound);
      msh_row = msh_evaluate_row(msh_bound,msh.nel_dir(2));
      msh_row.normal = msh_side.normal;                  % Import of structure normal.
      
      sp_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true, 'hessian', true, 'der3', true);
      % Normal and Tangent derivatives
      sp_normal = sp_normal_tang(sp_side_row, msh_row, 'gradient', true,'hessian', true, 'der3', true);
      
      sp_boundary = sp_side_row;
   end
   
   stiff_mat2 = stiff_mat2 + copp .* op_Bordo_Fracture(sp_boundary, sp_normal, msh_side);
 
end
stiff_mat = stiff_mat1 - stiff_mat2;
%   
u = zeros (space.ndof, 1);
% Apply Dirichlet boundary conditions (edge 3)
[u_drchlt3, drchlt_dofs3] = sp_drchlt_l2_proj (space, msh, h, 3);

iniz = ceil( size(drchlt_dofs3,1)/2 )-1;
fine = size(drchlt_dofs3,1);

drchlt_dofs3 = drchlt_dofs3(iniz : fine);
u_drchlt3 = u_drchlt3(iniz : fine);

% u_drchlt3=[]; 
% drchlt_dofs3=[];

% Apply Dirichlet boundary conditions (edges 1 2 4)
[u_drchlt2, drchlt_dofs2] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
% Aggiungo parte del bordo 3
drchlt_dofs = [drchlt_dofs2;drchlt_dofs3];
u_drchlt = [u_drchlt2; u_drchlt3];

u(drchlt_dofs) = u_drchlt;
%
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end