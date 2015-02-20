%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CALCOLA LE PARTI DI BORDO CHE NON VANNO MAI A ZERO                      %
% PER IL PROBLEMA ANTI-PLANE.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = op_Bordo_Fracture_grad (spu, spv, msh)

mat = spalloc (spv.ndof, spu.ndof, 1);

hessu = reshape (spu.shape_function_hessians, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
der3u = reshape (spu.shape_function_der3, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
funv  = reshape (spv.shape_functions, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);

ndirv = size (gradv, 2);
ndiru = size (hessu, 2);

for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
       mat_loc1 = zeros (spv.nsh(iel), spu.nsh(iel));
       mat_loc2 = zeros (spv.nsh(iel), spu.nsh(iel));
       for idof = 1:spv.nsh(iel)
           ishg = -reshape(gradv(:,:,:,idof,iel),spv.ncomp * ndirv, []);
           ishf = reshape(funv(:,:,:,idof,iel),spv.ncomp, []);
           for jdof = 1:spu.nsh(iel) 
               jshh = reshape(hessu(:,:,:,jdof,iel),spu.ncomp * ndiru, []);
               jshd = -reshape(der3u(:,:,:,jdof,iel),spu.ncomp * ndiru, []);
               % The cycle on the quadrature points is vectorized
               %for inode = 1:msh.nqn
               
               % Couple
               mat_loc1(idof, jdof) = mat_loc1(idof, jdof) + ...
                    sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ishg(2,:) .* jshh(1,:), 1).');   
               % Traction
               mat_loc2(idof, jdof) = mat_loc2(idof, jdof) + ... 
                   sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ishf(1,:) .* jshd(3,:), 1).');
                       
               %end
           end
       end
       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
           mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc1 + mat_loc2;
   else
       warning ('geopdes:jacdet_zero_at_quad_node',...
           'op_laplaceu_laplacev: singular map in element number %d', iel)
   end
end

end


% OLD FUNCTION!
% mat = spalloc(spu.ndof, spu.ndof, 1);
% 
% spv.connectivity = spu.connectivity;
% 
% for iel = 1:msh.nel
%     if (all (msh.jacdet(:,iel)))
%          mat_loc = zeros ( spu.nsh(iel), spu.nsh(iel) );
%          
%          mat_loc1 = zeros ( spu.nsh(iel), spu.nsh(iel) );
%          mat_loc2 = zeros ( spu.nsh(iel), spu.nsh(iel) );
%          
%          for idof = 1:spu.nsh(iel)
%              ispv = reshape (spu.shape_functions(:,idof,iel),[], msh.nqn); % shape function v
%              isnv = spv.deriv_normal(:,:,idof,iel);                        % shape deriv normal
%              
%              keyboard
%          for jdof = 1:spu.nsh(iel)
%              
%              jsnnu =  spv.deriv2_normal(:,:,jdof,iel);      % shape deriv2 normal
%              jsnnnu = spv.deriv3_normal(:,:,jdof,iel);      % shape deriv3 normal
%              
% % The cycle on the quadrature points is vectorized
%             %for inode = 1:msh.nqn
%             % TRACTION PART
%             mat_loc1(idof, jdof) = mat_loc1(idof, jdof) + ...
%               sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*sum ( jsnnnu(3,:) .* ispv,1).');
%             
%             % COUPLE PART
%             mat_loc2(idof, jdof) = mat_loc2(idof, jdof) + ...
%               sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*sum ( jsnnu(3,:).*isnv(1,:),1).'); 
%             
%             mat_loc = mat_loc1 + mat_loc2;
%            %end
%          end
%          end
%        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
%          mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) +  mat_loc;
%      else
%        warning ('geopdes:jacdet_zero_at_quad_node', 'op_Bordo_Couple: singular map in element number %d', iel)
%     end
% end
% 
% end