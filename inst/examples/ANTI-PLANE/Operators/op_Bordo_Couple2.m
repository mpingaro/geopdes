%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcola il pezzo al bordo dei momenti quando ho piastra appoggiata.    %
%                                                                         %
% (laplace_u , grad_u * n) on the boundary domain!!!!                     %
%                                                                         %
% Questo perché voglio solo imporre momenti nulli                         %
%                                                                         %
% Calcola l'integrale di bordo:  ( laplace(u) - d^2 u/dn^2) (du/dn)       %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat = op_Bordo_Couple2 (sp, msh, coeff)

mat = spalloc(sp.ndof, sp.ndof, 1);

%% NORMAL AND TANGENT VECTOR

% Comp. along x
nx = repmat(msh.normal(1,:,:), sp.nsh_max, 1);
nx = reshape(nx, msh.nqn, sp.nsh_max, msh.nel);
% Comp. along y
ny = repmat(msh.normal(2,:,:), sp.nsh_max, 1);
ny = reshape(ny, msh.nqn, sp.nsh_max, msh.nel);


% NORMAL AND TANGENT DERIVATIVES
u_x = reshape(sp.shape_function_gradients(1,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
u_y = reshape(sp.shape_function_gradients(2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
laplace = sp.shape_function_hessians(1,1,:,:,:)+ sp.shape_function_hessians(2,2,:,:,:);

laplace_x = sp.shape_function_der3(1,1,:,:,:) + sp.shape_function_der3(2,1,:,:,:);
laplace_y = sp.shape_function_der3(2,2,:,:,:) + sp.shape_function_der3(1,2,:,:,:);


laplace = reshape(laplace, msh.nqn, sp.nsh_max, msh.nel);
laplace_x = reshape(laplace_x, msh.nqn, sp.nsh_max, msh.nel);
laplace_y = reshape(laplace_y, msh.nqn, sp.nsh_max, msh.nel);

grad_n = u_x.*nx + u_y.*ny;                    % Normal Derivative
laplace_n = laplace_x.*nx + laplace_y.*ny;     % Normal of Laplace

spv  = reshape (sp.shape_functions, msh.nqn, sp.nsh_max, msh.nel);

for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
         mat_loc = zeros ( sp.nsh(iel), sp.nsh(iel) );
         mat_loc1 = zeros ( sp.nsh(iel), sp.nsh(iel) );
         mat_loc2 = zeros ( sp.nsh(iel), sp.nsh(iel) );
         mat_loc3 = zeros ( sp.nsh(iel), sp.nsh(iel) );
         
         for idof = 1:sp.nsh(iel)
             ipv = spv(:,idof,iel);      % funzione test
             ihv = laplace(:,idof,iel);
                
         for jdof = 1:sp.nsh(iel)
             jgnu = grad_n(:,jdof,iel);
             jhnu = laplace_n(:,jdof,iel);
             jhu = laplace(:,jdof,iel);
              
            %
            mat_loc1(idof, jdof) = mat_loc1(idof, jdof) + ...
              sum (msh.jacdet(:, iel).*msh.quad_weights(:, iel).*sum (jgnu.*ipv,1).');
            
            %
            mat_loc2(idof, jdof) = mat_loc2(idof, jdof) + ...
              sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel).*sum (jhnu.*ipv,1).'); 
            
            %
            mat_loc2(idof, jdof) = mat_loc2(idof, jdof) + ...
              sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*sum (jhu.*ihv,1).');
          
            mat_loc =  - mat_loc1 + (0.1*0.1).*(mat_loc2 - mat_loc3);
         end
         end
       mat(sp.connectivity(:, iel), sp.connectivity(:, iel)) = ...
         mat(sp.connectivity(:, iel), sp.connectivity(:, iel)) +  mat_loc;
     else
       warning ('geopdes:jacdet_zero_at_quad_node', 'op_Bordo_Couple: singular map in element number %d', iel)
    end
end

end