function sp_side = sp_eval_boundary_col (space, msh, normal, varargin)

%iside = msh_side.side_number;
%ind = mod (floor ((iside+1)/2), 2) + 1;

%bnodes = reshape (squeeze (msh_side.quad_nodes(ind,:,:)), msh_side.nqn, []);

%sp_side = sp_bspline_1d_param (sp.knots{ind}, sp.degree(ind), bnodes, ...
%    'gradient', true, 'hessian', true);

% switch (iside)
%     case 1
%         weights = sp.weights(1,:);
%    case 2
%         weights = sp.weights(end,:);
%    case 3
%         weights = sp.weights(:,1);
%    case 4
%         weights = sp.weights(:,end);
%  end
% 
%   sp_side = bsp_2_nrb_1d__ (sp_side, msh_side, weights);

%  sp_side.dofs = sp.boundary(iside).dofs;
  

%% VARARGIN ----> Valutazione == Valore,  Gradinets, Hessian.   
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


sp_side = sp_eval_boundary_col_param (space, msh, normal, varargin{:});
%%
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
      bu = reshape (sp_side.shape_function_gradients(1,:,ii,:), msh.nqn, msh.nel);
      bv = reshape (sp_side.shape_function_gradients(2,:,ii,:), msh.nqn, msh.nel);
      buu = reshape (sp_side.shape_function_hessians(1,1,:,ii,:), msh.nqn, msh.nel);
      buv = reshape (sp_side.shape_function_hessians(1,2,:,ii,:), msh.nqn, msh.nel);
      bvv = reshape (sp_side.shape_function_hessians(2,2,:,ii,:), msh.nqn, msh.nel);
        
      [bxx, bxy, byy] = der2_basisfun_phys__ (xu(:), xv(:), yu(:),...
          yv(:), xuu(:), xuv(:), xvv(:), yuu(:), yuv(:), yvv(:),...
          buu(:), buv(:), bvv(:),bu(:),bv(:));
      
      sh = size (sp.shape_function_hessians(1,1,:,ii,:));
      sp_side.shape_function_hessians(1,1,:,ii,:) = reshape (bxx, sh);
      sp_side.shape_function_hessians(1,2,:,ii,:) = reshape (bxy, sh);
      sp_side.shape_function_hessians(2,1,:,ii,:) = reshape (bxy, sh);
      sp_side.shape_function_hessians(2,2,:,ii,:) = reshape (byy, sh);
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
      bu = reshape (sp_side.shape_function_gradients(1,:,ii,:), msh.nqn, msh.nel);
      bv = reshape (sp_side.shape_function_gradients(2,:,ii,:), msh.nqn, msh.nel);
      buu = reshape (sp_side.shape_function_hessians(1,1,:,ii,:), msh.nqn, msh.nel);
      buv = reshape (sp_side.shape_function_hessians(1,2,:,ii,:), msh.nqn, msh.nel);
      bvv = reshape (sp_side.shape_function_hessians(2,2,:,ii,:), msh.nqn, msh.nel);
      buuu = reshape (sp_side.shape_function_der3(1,1,:,ii,:), msh.nqn, msh.nel);
      buuv = reshape (sp_side.shape_function_der3(1,2,:,ii,:), msh.nqn, msh.nel);
      buvv = reshape (sp_side.shape_function_der3(2,1,:,ii,:), msh.nqn, msh.nel);
      bvvv = reshape (sp_side.shape_function_der3(2,2,:,ii,:), msh.nqn, msh.nel);
        
      [bxxx, bxxy, bxyy, byyy] = der3_basisfun_phys__ (xu, xv, yu,...
          yv, xuu, xuv, xvv, yuu, yuv, yvv, xuuu, xuuv, xuvv, xvvv,...
          yuuu, yuuv, yuvv, yvvv, buuu, buuv, buvv, bvvv, buu, buv, bvv, bu, bv);
      
      sh = size (sp.shape_function_der3(1,:,ii,:));
      sp_side.shape_function_der3(1,1,:,ii,:) = reshape (bxxx, sh);
      sp_side.shape_function_der3(1,2,:,ii,:) = reshape (bxxy, sh);
      sp_side.shape_function_der3(2,1,:,ii,:) = reshape (bxyy, sh);
      sp_side.shape_function_der3(2,2,:,ii,:) = reshape (byyy, sh);
  end
  clear bu bv buu buv bvv buuu buuv buvv bvvv
  clear xu xv yu yv xuu xuv xvv yuu yuv yvv xuuu xuuv xuvv xvvv yuuu yuuv yuvv yvvv ...
      bxxx bxxy bxyy byyy
end

if (gradient)
  JinvT = geopdes_invT__ (msh.geo_map_jac);
  JinvT = reshape (JinvT, [2, 2, msh.nqn, msh.nel]);
  sp_side.shape_function_gradients = geopdes_prod__ (JinvT, sp_side.shape_function_gradients);
end

end