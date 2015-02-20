% OP_FRACTURE_ANTIPLANEGRAD_TP: assemble the stiffness matrix 
% A = [a(i,j)], a(i,j) = (ls^2 laplace u_j, laplace v_i) + (grad u_j,grad v_i) , exploiting the tensor product structure.
%
%   mat = op_fracture_antiPlane_tp (spu, spv, msh, epsilon);
%   [rows, cols, values] = op_fracture_antiPlane_tp (spu, spv, msh, epsilon);
%
% INPUT:
%
%   spu:     class representing the space of trial functions (see sp_bspline_2d)
%   spv:     class representing the space of test functions (see sp_bspline_2d)
%   msh:     class defining the domain partition and the quadrature rule (see msh_2d)
%   epsilon: function handle to compute the diffusion coefficient
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
function varargout = op_fracture_antiPlaneGRAD_tp (space1, space2, msh, coeff1, coeff2)

  A1 = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);
  A2 = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);
  
  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true, 'hessian', true, 'der3', false);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true, 'hessian', true, 'der3', false);
    
    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A1 = A1 + op_gradu_gradv (sp1_col, sp2_col, msh_col, coeff1(x{:}));
    A2 = A2 + op_gradgradu_gradgradv (sp1_col, sp2_col , msh_col, coeff2(x{:}));
  end
  
  A = A1 + A2; 

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end