function varargout = op_Bordo_Couple_tp (space1, space2, msh)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  ndim = numel (msh.qn);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', true, 'gradient', true, 'hessian', true, 'der3', false);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', true, 'gradient', true, 'hessian', true, 'der3', false);
    
    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_Bordo_Couple (sp1_col,sp2_col, msh_col);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end