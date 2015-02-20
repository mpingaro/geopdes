% MSH_EVALUATE_ROW: evaluate the parameterization in one row of the mesh.
%
%     msh_row = msh_evaluate_row (msh, rownum)
%
% INPUTS:
%
%    msh:    mesh object (see msh_2d)
%    rownum: number of the "row", i.e., the element in the second parametric direction.
%
% OUTPUT:
%
%     msh_col: structure containing the quadrature rule in one column of the physical domain, which contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     rownum        (scalar)                number of the row
%     nel           (scalar)                number of elements in the row
%     elem_list     (nel vector)            indices of the elements in the row
%     nel_dir       (1 x 2 vector)          number of elements in each parametric direction for the entire mesh
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 2 vector)          number of quadrature nodes per element in each parametric direction
%     quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric space
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%     geo_map_der2  (2 x 2 x 2 x nqn x nel) Second order derivatives of the map evaluated at the quadrature nodes (optional)
%     jacdet        (nqn x nel)             determinant of the Jacobian evaluated in the quadrature points

function msh_row = msh_evaluate_row (msh, rownum)

  msh_row.rownum = rownum;
  msh_row.elem_list = msh.nel_dir(1)*(rownum-1)*ones(1,msh.nel_dir(1))+(1:msh.nel_dir(1));

  msh_row.nel_dir = msh.nel_dir;
  msh_row.nel  = msh.nel_dir(1);

  msh_row.nqn     = msh.nqn;
  msh_row.nqn_dir = msh.nqn_dir;

  qnu = msh.qn{1}; qnv = msh.qn{2}(:,rownum);
  
  if (isempty (msh.quad_nodes))
    quad_nodes_u = qnu;
    quad_nodes_u = repmat  (quad_nodes_u, [msh.nqn_dir(2), 1, 1]);
    quad_nodes_u = reshape (quad_nodes_u, [], msh.nel_dir(1));
    
    quad_nodes_v = reshape (qnv, msh.nqn_dir(2), 1, 1);
    quad_nodes_v = repmat  (quad_nodes_v, [1, msh.nqn_dir(1), msh.nel_dir(1)]);
    quad_nodes_v = reshape (quad_nodes_v, [], msh.nel_dir(1));
    
    msh_row.quad_nodes(1, :, :) = quad_nodes_u;
    msh_row.quad_nodes(2, :, :) = quad_nodes_v;

    clear quad_nodes_u quad_nodes_v
  else
    msh_row.quad_nodes = msh.quad_nodes(:,:,msh_row.elem_list);
  end


  if (~isempty (msh.qw))
    if (isempty (msh.quad_weights))
      qwu = msh.qw{1};  qwv = msh.qw{2}(:,rownum);
      
      quad_weights_u = qwu;
      quad_weights_u = repmat  (quad_weights_u, [msh.nqn_dir(2), 1, 1]);
      quad_weights_u = reshape (quad_weights_u, [], msh.nel_dir(1));
      
      quad_weights_v = reshape (qwv, msh.nqn_dir(2), 1, 1);
      quad_weights_v = repmat  (quad_weights_v, [1, msh.nqn_dir(1), msh.nel_dir(1)]);
      quad_weights_v = reshape (quad_weights_v, [], msh.nel_dir(1));
      
      msh_row.quad_weights = quad_weights_u .* quad_weights_v;

      clear quad_weights_u quad_weights_v
    else
      msh_row.quad_weights = msh.quad_weights(:,msh_row.elem_list);
    end
  end

  if (isempty (msh.geo_map))
    F = feval (msh.map, {qnu(:)', qnv(:)'});
    msh_row.geo_map = reshape (F, [2, msh.nqn, msh.nel_dir(1)]);
  else
    msh_row.geo_map = msh.geo_map(:,:,msh_row.elem_list);
  end
  
  if (isempty (msh.geo_map_jac))
    jac = feval (msh.map_der, {qnu(:)', qnv(:)'});
    msh_row.geo_map_jac = reshape (jac, 2, 2, msh.nqn, msh.nel_dir(1));
  else
    msh_row.geo_map_jac = msh.geo_map_jac(:,:,:,msh_row.elem_list);
  end
   
  if (isempty (msh.jacdet))
    msh_row.jacdet = abs (geopdes_det__ (msh_row.geo_map_jac));
    msh_row.jacdet = reshape (msh_row.jacdet, [msh.nqn, msh.nel_dir(1)]);
  else
    msh_row.jacdet = msh.jacdet(:,msh_row.elem_list);
  end

  if (msh.der2)
    if (isempty (msh.geo_map_der2))
      msh_row.geo_map_der2 = reshape (feval (msh.map_der2, {qnu(:)', qnv(:)'}), ...
                                         2, 2, 2, msh_row.nqn, msh_row.nel);
    else
      msh_row.geo_map_der2 = msh.geo_map_der2(:,:,:,:,msh_row.elem_list);
    end
  end
  
  if (msh.der3)
    if (isempty (msh.geo_map_der3))
      msh_row.geo_map_der3 = reshape (feval (msh.map_der3, {qnu(:)', qnv(:)'}), ...
                                         2, 2, 2, msh_row.nqn, msh_row.nel);
    else
      msh_row.geo_map_der3 = msh.geo_map_der3(:,:,:,:,msh_row.elem_list);
      
    end
  end
  
end