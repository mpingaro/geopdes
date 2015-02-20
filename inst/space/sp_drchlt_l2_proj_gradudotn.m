% SP_DRCHLT_L2_PROJ_GRADUDOTN: assign the degrees of freedom of Dirichlet boundaries (rotations du/dn) through an L2 projection.
%
%   [u, dofs_ad] = sp_drchlt_l2_proj_gradudotn (sp, msh, ud, g, sides)
%
% INPUT:
%
%  sp:       object defining the space of discrete functions (see sp_bspline_2d)
%  msh:      object defining the domain partition and the quadrature rule (see msh_2d)
%  geometry: 
%  ud:       vector contains the values of displacement Dirichlet condition
%  g:        function handle to compute the Dirichlet condition (rotations)
%  sides:    boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:       assigned value to the degrees of freedom
%  dofs_ad: global numbering of the corresponding basis functions
%
% Copyright (C) 2010, Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, Rafael Vazquez
% Copyright (C) 2015, Rafael Vazquez, Marco Pingaro
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [u, dofs_ad] = sp_drchlt_l2_proj_gradudotn (sp, msh, geometry, ud, g, sides)

  rhs  = zeros (sp.ndof, 1);
  M = spalloc (sp.ndof, sp.ndof, 3*sp.ndof);

  dofs = [];
  dofs_ad = [];
  nent = 0;
  for iside = sides
      nent = nent + msh.boundary(iside).nel * sp.boundary(iside).nsh_max^2;
      dofs_ad = union (dofs_ad, sp.boundary(iside).adjacent_dofs);
      dofs    = union (dofs, sp.boundary(iside).dofs);
  end
  dofs_ad = setdiff(dofs_ad,dofs);
  
  for iside = sides
      msh_bnd = msh_eval_boundary_side (msh, iside);
      x = squeeze (msh_bnd.geo_map(1,:,:));
      y = squeeze (msh_bnd.geo_map(2,:,:));
      gval = reshape (g (x, y, iside), msh_bnd.nqn, msh_bnd.nel);
      
      if (iside == 1)
        breaks = {msh.breaks{1}(1:2) msh.breaks{2}};
        qn = {msh.breaks{1}(1), msh.qn{2}};
        qw = {1, msh.qw{2}};
        msh_aux = msh_2d (breaks, qn, qw, geometry, 'der2', true);
        space_side = sp_nurbs_2d (geometry.nurbs, msh_aux);
        
        msh_col = msh_evaluate_col(msh_aux,1);
        msh_col.normal = msh_bnd.normal;     
        
        space_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true);
        space_normal = sp_normal_tang(space_side_col, msh_col, 'gradient', true);
        
        space_side.boundary(1).dofs = space_side.boundary(1).adjacent_dofs;
        sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true);
        sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true);
        
        rhs = rhs + op_f_gradu_n (sp_normal, msh_col, gval);
        rhs_ad = op_u_norm_v(sp_normal, space_normal, msh_col, ones (msh_col.nqn, msh_col.nel))*ud;
        rhs = rhs - rhs_ad;
        
        M = M + op_u_norm_v(sp_normal, sp_normal, msh_col, ones (msh_col.nqn, msh_col.nel));
        
      elseif (iside == 2)
        breaks = {msh.breaks{1}(end-1:end) msh.breaks{2}};
        qn = {msh.breaks{1}(end), msh.qn{2}};
        qw = {1, msh.qw{2}};
        msh_aux = msh_2d (breaks, qn, qw, geometry, 'der2', true);
        space_side = sp_nurbs_2d (geometry.nurbs, msh_aux);
        
        msh_col = msh_evaluate_col(msh_aux,msh_aux.nel_dir(1));
        msh_col.normal = msh_bnd.normal;
        
        space_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true);
        space_normal = sp_normal_tang(space_side_col, msh_col, 'gradient', true);
        
        space_side.boundary(2).dofs = space_side.boundary(2).adjacent_dofs;
        sp_side_col = sp_evaluate_col(space_side,msh_col,'gradient', true); 
        sp_normal = sp_normal_tang(sp_side_col, msh_col, 'gradient', true);
        
        rhs = rhs + op_f_gradu_n (sp_normal, msh_col, gval);
        rhs_ad = op_u_norm_v(sp_normal, space_normal, msh_col, ones (msh_col.nqn, msh_col.nel))*ud;
        rhs = rhs - rhs_ad; 
        
        M = M + op_u_norm_v(sp_normal, sp_normal, msh_col, ones (msh_col.nqn, msh_col.nel));
        
      elseif (iside == 3)
        breaks = {msh.breaks{1} msh.breaks{2}(1:2)};
        qn = {msh.qn{1}, msh.breaks{2}(1)};
        qw = {msh.qw{1}, 1};
        msh_aux = msh_2d (breaks, qn, qw, geometry, 'der2', true);
        space_side = sp_nurbs_2d (geometry.nurbs, msh_aux);
        
        msh_row = msh_evaluate_row(msh_aux,1);
        msh_row.normal = msh_bnd.normal;
        
        space_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true);
        space_normal = sp_normal_tang(space_side_row, msh_row, 'gradient', true);
        
        space_side.boundary(3).dofs = space_side.boundary(3).adjacent_dofs;
        sp_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true);
        sp_normal = sp_normal_tang(sp_side_row, msh_row, 'gradient', true);
        
        rhs = rhs + op_f_gradu_n (sp_normal, msh_row, gval);
        rhs_ad = op_u_norm_v(sp_normal, space_normal, msh_row, ones (msh_row.nqn, msh_row.nel))*ud;
        rhs = rhs - rhs_ad;
        
        M = M + op_u_norm_v(sp_normal, sp_normal, msh_row, ones (msh_row.nqn, msh_row.nel));
        
      elseif (iside == 4)
        breaks = {msh.breaks{1} msh.breaks{2}(end-1:end)};
        qn = {msh.qn{1}, msh.breaks{2}(end)};
        qw = {msh.qw{1}, 1};
        msh_aux = msh_2d (breaks, qn, qw, geometry, 'der2', true);
        space_side = sp_nurbs_2d (geometry.nurbs, msh_aux);
        
        msh_row = msh_evaluate_row(msh_aux,msh_aux.nel_dir(2));
        msh_row.normal = msh_bnd.normal;
        
        space_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true);
        space_normal = sp_normal_tang(space_side_row, msh_row, 'gradient', true);
        
        space_side.boundary(4).dofs = space_side.boundary(4).adjacent_dofs;
        sp_side_row = sp_evaluate_row(space_side,msh_row,'gradient', true);
        sp_normal = sp_normal_tang(sp_side_row, msh_row, 'gradient', true);

        rhs = rhs + op_f_gradu_n (sp_normal, msh_bnd, gval);
        rhs_ad = op_u_norm_v(sp_normal, space_normal, msh_row, ones (msh_row.nqn, msh_row.nel))*ud;
        rhs = rhs - rhs_ad;
        
        M = M + op_u_norm_v(sp_normal, sp_normal, msh_row, ones (msh_row.nqn, msh_row.nel));
        
      end      
  end
  
  u = M(dofs_ad, dofs_ad) \ rhs(dofs_ad, 1);
end