function sp_side = sp_eval_boundary_col_param (space, normal, varargin)

value = true;
gradient = false;
hessian = false;
der3 = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col_param: options must be passed in the [option, value] format');
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
      error ('sp_evaluate_col_param: unknown option %s', varargin {ii});
    end
  end
end

% NORMAL AND TANGENT VECTOR ( ONLY PLANE CASES!!! )
norm_x = normal(1,:,:);
norm_y = normal(2,:,:);
tang_x = -normal(2,:,:);
tang_y = normal(1,:,:);
  
if (gradient)
  % NORMAL AND TANGENT DERIVATIVES
  sp_side.deriv_normal(1,:,:,:) = spu.shape_function_gradients(:,:,:) .* norm_x +...
      spu.shape_function_gradients(:,:,:) .* norm_y; % Normal Derivative
  sp_side.deriv_normal(2,:,:,:) = spu.shape_function_gradients(:,:,:) .* norm_y +...
      spu.shape_function_gradients(:,:,:) .* norm_x; % Tangent Derivative
end
                           
if (hessian)
      % SECOND NORMAL AND TANGENT DERIVATIVES
      shape_function_hessians = space.shape_function_hessians;
      
      sp_side.deriv2_normal(1,:,:,:) = shape_function_hessians(1,1,:,:,:) .* norm_x.^2 +...
          2*shape_function_hessians(1,2,:,:,:) .* norm_x .* norm_y + ...
          shape_function_hessians(2,2,:,:,:) .* norm_y.^2;
      
      sp_side.deriv2_normal(2,:,:,:) = shape_function_hessians(1,1,:,:,:) .* norm_x .* tang_x +...
          shape_function_hessians(1,2,:,:,:) .* norm_x .* tang_y + ...
          shape_function_hessians(2,1,:,:,:) .* norm_y .* tang_x + ...
          shape_function_hessians(2,2,:,:,:) .* norm_y .* tang_y;
      
      sp_side.deriv2_normal(3,:,:,:) = shape_function_hessians(1,1,:,:,:) .* tang_x.^2 +...
          2*shape_function_hessians(1,2,:,:,:) .* tang_x .* tang_y + ...
          shape_function_hessians(2,2,:,:,:) .* tang_y.^2;   
end
  
if(der3)
      shape_function_der3 = space.shape_function_der3;
      
      % THIRD NORMAL AND TANGENT DERIVATIVES
      % R_nnn
      sp_side.deriv3_normal(4,:,:,:) = shape_function_der3(1,1,:,:,:) .* norm_x.^3 +...
          3*shape_function_der3(1,2,:,:,:) .* norm_x .* norm_x .* norm_y +...
          3*shape_function_der3(2,1,:,:,:) .* norm_x .* norm_y .* norm_y +...
          shape_function_der3(2,2,:,:,:) .* norm_y.^3;
      
      % R_nnt & R_ntn % R_tnn
      sp_side.deriv3_normal(2,:,:,:) = shape_function_der3(1,1,:,:,:) .* norm_x .* norm_x .* tang_x +...
          3*shape_function_der3(1,2,:,:,:) .* (norm_x .* norm_x .* tang_y + norm_x .* norm_y .* tang_x + norm_y .* norm_x .* tang_x) +...
          3*shape_function_der3(2,1,:,:,:) .* (norm_x .* norm_y .* tang_y + norm_y .* norm_x .* tang_y + norm_y .* norm_y .* tang_x) +...
          shape_function_der3(2,2,:,:,:) .* norm_y .* norm_y .* tang_y;

      % R_ntt & R_ttn & R_tnt
      sp_side.deriv3_normal(3,:,:,:) = shape_function_der3(1,1,:,:,:) .* norm_x .* tang_x .* tang_x +...
          3*shape_function_der3(1,2,:,:,:) .* (norm_x .* tang_x .* tang_y + norm_x .* tang_y .* tang_x + norm_y .* tang_x .* tang_x) +...
          3*shape_function_der3(2,1,:,:,:) .* (norm_x .* tang_y .* tang_y + norm_y .* tang_x .* tang_y + norm_y .* tang_y .* tang_x) +...
          shape_function_der3(2,2,:,:,:) .* norm_y .* tang_y .* tang_y;
      
      % R_ttt
      sp_side.deriv3_normal(4,:,:,:) = shape_function_der3(1,1,:,:,:) .* tang_x.^3 +...
          3*shape_function_der3(1,2,:,:,:) .* tang_x .* tang_x .* tang_y +...
          3*shape_function_der3(2,1,:,:,:) .* tang_x .* tang_y .* tang_y +...
          shape_function_der3(2,2,:,:,:) .* tang_y.^3;
end
   
end