function sp_side = sp_normal_tang(space, msh, varargin)

gradient = false;
hessian = false;
der3 = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col_param: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
        hessian = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'der3'))
        der3 = varargin {ii+1};
    else
      error ('sp_evaluate_col_param: unknown option %s', varargin {ii});
    end
  end
end

%% NORMAL AND TANGENT VECTOR

% Comp. along x
normal_x = repmat(msh.normal(1,:,:), space.nsh_max, 1);
normal_x = reshape(normal_x, msh.nqn, space.nsh_max, msh.nel);

% Comp. along y
normal_y = repmat(msh.normal(2,:,:), space.nsh_max, 1);
normal_y = reshape(normal_y, msh.nqn, space.nsh_max, msh.nel);
% Tangential vector
tang_x = -normal_y;
tang_y = normal_x;
%%
if (gradient)
  % NORMAL AND TANGENT DERIVATIVES
  grad_x = reshape(space.shape_function_gradients(1,:,:,:), msh.nqn, space.nsh_max, msh.nel);
  grad_y = reshape(space.shape_function_gradients(2,:,:,:), msh.nqn, space.nsh_max, msh.nel);
  
  sp_side.deriv_normal(1,:,:,:) = grad_x .* normal_x + grad_y .* normal_y; % Normal Derivative
  sp_side.deriv_normal(2,:,:,:) = grad_x .* tang_x + grad_y .* tang_y;     % Tangent Derivative
end
                           
if (hessian)
    % SECOND NORMAL AND TANGENT DERIVATIVES
    hessian_xx = reshape(space.shape_function_hessians(1,1,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    hessian_xy = reshape(space.shape_function_hessians(1,2,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    hessian_yy = reshape(space.shape_function_hessians(2,2,:,:,:), msh.nqn, space.nsh_max, msh.nel);
     
    % R_nn
    sp_side.deriv2_normal(1,:,:,:) = hessian_xx.*normal_x.^2 +...
        2.*hessian_xy.*normal_x.*normal_y + hessian_yy.*normal_y.^2;
    % R_nt & R_tn
    sp_side.deriv2_normal(2,:,:,:) = hessian_xx .* normal_x .* tang_x +...
        hessian_xy.*(normal_x.*tang_y + normal_y.*tang_x) + ...
        hessian_yy.*normal_y.*tang_y;
    % R_tt
    sp_side.deriv2_normal(3,:,:,:) = hessian_xx.*tang_x.^2 +...
        2.*hessian_xy.*tang_x.*tang_y + hessian_yy.*tang_y.^2;  
end
  
if(der3)
    der3_xxx = reshape(space.shape_function_der3(1,1,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    der3_xxy = reshape(space.shape_function_der3(1,2,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    der3_xyy = reshape(space.shape_function_der3(2,1,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    der3_yyy = reshape(space.shape_function_der3(2,2,:,:,:), msh.nqn, space.nsh_max, msh.nel);
    
    % THIRD NORMAL AND TANGENT DERIVATIVES
    % R_nnn
    sp_side.deriv3_normal(1,:,:,:) = der3_xxx.*normal_x.^3 +...
        3.*der3_xxy.*normal_x.^2.* normal_y + 3.*der3_xyy.*normal_x.*normal_y.^2 +...
        der3_yyy.*normal_y.^3;
    
    % R_nnt & R_ntn % R_tnn
    sp_side.deriv3_normal(2,:,:,:) = der3_xxx .* normal_x .^2.* tang_x +...
        der3_xxy .* (normal_x.^2.*tang_y + 2.*normal_x.*normal_y.*tang_x) +...
        der3_xyy .* (2.*normal_x.*normal_y.*tang_y + normal_y.^2.*tang_x) +...
        der3_yyy.*normal_y.^2.*tang_y;
    
    % R_ntt & R_ttn & R_tnt
    sp_side.deriv3_normal(3,:,:,:) = der3_xxx.*tang_x.^2.*normal_x +...
        der3_xxy.*(tang_x.^2.*normal_y+ 2.*tang_x.*tang_y.*normal_x) +...
        der3_xyy.*(2.*tang_x.*tang_y.*normal_y + tang_y.^2.*normal_x) +...
        der3_yyy.*tang_y.^2.*normal_y;
    
    % R_ttt
    sp_side.deriv3_normal(4,:,:,:) = der3_xxx.*tang_x.^3 +...
        3.*der3_xxy.*tang_x.^2.*tang_y + 3.*der3_xyy .* tang_x .* tang_y.^2 +...
        der3_yyy.*tang_y.^3;
end

% Copy structure of space
sp_side.nsh_max = space.nsh_max;
sp_side.nsh = space.nsh;
sp_side.ndof = space.ndof;
sp_side.ndof_dir = space.ndof_dir;
sp_side.connectivity = space.connectivity;
sp_side.ncomp = space.ncomp;

end