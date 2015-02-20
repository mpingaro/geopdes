% SP_EVALUATE_ROW: compute the basis functions in one column of the mesh.
%
%     sp = sp_evaluate_row (space, msh_row, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:   object defining the space of discrete functions (see sp_nurbs_2d)
%    msh_row: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_2d/msh_evaluate_row)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%            der3       |      false      |  compute shape_function_der3
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                            DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (1 x 2 vector)                         degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh_row.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh_row.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_row.nqn x nsh_max x msh_row.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%                (2 x msh_row.nqn x nsh_max x msh_row.nel)  basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%             (2 x 2 x msh_row.nqn x nsh_max x msh_row.nel) basis function hessians evaluated at each quadrature node in each element
%    shape_function_der3
%             (2 x 2 x msh_row.nqn x nsh_max x msh_row.nel)  third derivatives of basis functions evaluated at each quadrature node in each element

function sp = sp_evaluate_row (space, msh, varargin)

value = true;
gradient = false;
hessian = false;
der3 = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
        hessian = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'der3'))
        der3 = varargin {ii+1};
    else
      error ('sp_evaluate_col: unknown option %s', varargin {ii});
    end
  end
end

sp = sp_evaluate_row_param (space, msh, varargin{:});

if (hessian && isfield (msh, 'geo_map_der2'))
  xu = reshape (msh.geo_map_jac(1,1,:,:), msh.nqn, msh.nel);
  xv = reshape (msh.geo_map_jac(1,2,:,:), msh.nqn, msh.nel);
  yu = reshape (msh.geo_map_jac(2,1,:,:), msh.nqn, msh.nel);
  yv = reshape (msh.geo_map_jac(2,2,:,:), msh.nqn, msh.nel);

  xuu = reshape (msh.geo_map_der2(1, 1, 1, :, :), [], msh.nel);
  yuu = reshape (msh.geo_map_der2(2, 1, 1, :, :), [], msh.nel);
  xuv = reshape (msh.geo_map_der2(1, 1, 2, :, :), [], msh.nel);
  yuv = reshape (msh.geo_map_der2(2, 2, 1, :, :), [], msh.nel);
  xvv = reshape (msh.geo_map_der2(1, 2, 2, :, :), [], msh.nel);
  yvv = reshape (msh.geo_map_der2(2, 2, 2, :, :), [], msh.nel);

    for ii=1:sp.nsh_max   
        bu = reshape (sp.shape_function_gradients(1,:,ii,:), msh.nqn, msh.nel);
        bv = reshape (sp.shape_function_gradients(2,:,ii,:), msh.nqn, msh.nel);
        buu = reshape (sp.shape_function_hessians(1,1,:,ii,:), msh.nqn, msh.nel);
        buv = reshape (sp.shape_function_hessians(1,2,:,ii,:), msh.nqn, msh.nel);
        bvv = reshape (sp.shape_function_hessians(2,2,:,ii,:), msh.nqn, msh.nel);
        
        [bxx, bxy, byy] = der2_basisfun_phys__ (xu(:), xv(:), yu(:),...
            yv(:), xuu(:), xuv(:), xvv(:), yuu(:), yuv(:), yvv(:),...
            buu(:), buv(:), bvv(:),bu(:),bv(:));

        sh = size (sp.shape_function_hessians(1,1,:,ii,:));
        sp.shape_function_hessians(1,1,:,ii,:) = reshape (bxx, sh);
        sp.shape_function_hessians(1,2,:,ii,:) = reshape (bxy, sh);
        sp.shape_function_hessians(2,1,:,ii,:) = reshape (bxy, sh);
        sp.shape_function_hessians(2,2,:,ii,:) = reshape (byy, sh);
    end
    clear bu bv buu buv bvv 
    clear xu xv yu yv xuu xuv xvv yuu yuv yvv bxx bxy byy

end


if (der3 && isfield (msh, 'geo_map_der3'))
  xu = reshape (msh.geo_map_jac(1,1,:,:), msh.nqn, msh.nel);
  xv = reshape (msh.geo_map_jac(1,2,:,:), msh.nqn, msh.nel);
  yu = reshape (msh.geo_map_jac(2,1,:,:), msh.nqn, msh.nel);
  yv = reshape (msh.geo_map_jac(2,2,:,:), msh.nqn, msh.nel);

  xuu = reshape (msh.geo_map_der2(1, 1, 1, :, :), [], msh.nel);
  yuu = reshape (msh.geo_map_der2(2, 1, 1, :, :), [], msh.nel);
  xuv = reshape (msh.geo_map_der2(1, 1, 2, :, :), [], msh.nel);
  yuv = reshape (msh.geo_map_der2(2, 2, 1, :, :), [], msh.nel);
  xvv = reshape (msh.geo_map_der2(1, 2, 2, :, :), [], msh.nel);
  yvv = reshape (msh.geo_map_der2(2, 2, 2, :, :), [], msh.nel);
  
  xuuu = reshape (msh.geo_map_der3(1, 1, 1, :, :), [], msh.nel);
  yuuu = reshape (msh.geo_map_der3(2, 1, 1, :, :), [], msh.nel);
  xuuv = reshape (msh.geo_map_der3(1, 1, 2, :, :), [], msh.nel);
  yuuv = reshape (msh.geo_map_der3(2, 1, 2, :, :), [], msh.nel);
  xuvv = reshape (msh.geo_map_der3(1, 2, 1, :, :), [], msh.nel);
  yuvv = reshape (msh.geo_map_der3(2, 2, 1, :, :), [], msh.nel);
  xvvv = reshape (msh.geo_map_der3(1, 2, 2, :, :), [], msh.nel);
  yvvv = reshape (msh.geo_map_der3(2, 2, 2, :, :), [], msh.nel);
  

  for ii=1:sp.nsh_max   
        bu = reshape (sp.shape_function_gradients(1,:,ii,:), msh.nqn, msh.nel);
        bv = reshape (sp.shape_function_gradients(2,:,ii,:), msh.nqn, msh.nel);
        buu = reshape (sp.shape_function_hessians(1,1,:,ii,:), msh.nqn, msh.nel);
        buv = reshape (sp.shape_function_hessians(1,2,:,ii,:), msh.nqn, msh.nel);
        bvv = reshape (sp.shape_function_hessians(2,2,:,ii,:), msh.nqn, msh.nel);
        buuu = reshape (sp.shape_function_der3(1,1,:,ii,:), msh.nqn, msh.nel);
        buuv = reshape (sp.shape_function_der3(1,2,:,ii,:), msh.nqn, msh.nel);
        buvv = reshape (sp.shape_function_der3(2,1,:,ii,:), msh.nqn, msh.nel);
        bvvv = reshape (sp.shape_function_der3(2,2,:,ii,:), msh.nqn, msh.nel);
        
        [bxxx, bxxy, bxyy, byyy] = der3_basisfun_phys__ (xu, xv, yu,...
            yv, xuu, xuv, xvv, yuu, yuv, yvv, xuuu, xuuv, xuvv, xvvv,...
            yuuu, yuuv, yuvv, yvvv, buuu, buuv, buvv, bvvv, buu, buv, bvv, bu, bv);

        sh = size (sp.shape_function_der3(1,1,:,ii,:));
        sp.shape_function_der3(1,1,:,ii,:) = reshape (bxxx, sh);
        sp.shape_function_der3(1,2,:,ii,:) = reshape (bxxy, sh);
        sp.shape_function_der3(2,1,:,ii,:) = reshape (bxyy, sh);
        sp.shape_function_der3(2,2,:,ii,:) = reshape (byyy, sh);
  end
    clear bu bv buu buv bvv buuu buuv buvv bvvv 
    clear xu xv yu yv xuu xuv xvv yuu yuv yvv xuuu xuuv xuvv xvvv yuuu yuuv yuvv yvvv ...
          bxxx bxxy bxyy byyy
end

if (gradient)
  JinvT = geopdes_invT__ (msh.geo_map_jac);
  JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
  sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
end

end