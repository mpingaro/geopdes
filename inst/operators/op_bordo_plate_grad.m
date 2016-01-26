function mat = op_bordo_plate_grad (spu, spv, msh, pt)

mat = spalloc (spv.ndof, spu.ndof, 1);
for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
       mat_loc_mn = zeros (spv.nsh(iel), spu.nsh(iel));
       mat_loc_mnt = zeros (spv.nsh(iel), spu.nsh(iel));
       mat_loc_v_1 = zeros (spv.nsh(iel), spu.nsh(iel));
       mat_loc_v_2 = zeros (spv.nsh(iel), spu.nsh(iel));
       for idof = 1:spv.nsh(iel)
           ishh = spv.deriv2_normal(:,:,idof,iel);
           % ish3 = spv.deriv3_normal(:,:,idof, iel);
           for jdof = 1:spu.nsh(iel) 
               jshg = spv.deriv_normal(:,:,jdof,iel);
               % jshp = spu.shape_functions(:,jdof, iel);
               
               % Couple Mn  
               mat_loc_mn(idof, jdof) = mat_loc_mn(idof, jdof) + ...
                    sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ishh(3,:) .* jshg(1,:), 1).'); 
               % Couple Mnt
               mat_loc_mnt(idof, jdof) = mat_loc_mnt(idof, jdof) + ...
                    sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ishh(2,:) .* jshg(2,:), 1).');
%                Shear -> Third derivatives
%                mat_loc_v_1(idof, jdof) = mat_loc_v_1(idof, jdof) + ... 
%                    sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ish3(1,:) .* jshp', 1).');
%                mat_loc_v_2(idof, jdof) = mat_loc_v_2(idof, jdof) + ... 
%                    sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* sum( ish3(3,:) .* jshp', 1).');
                       
           end
       end
       mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
           mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) +...
           pt.*mat_loc_mn - pt.*mat_loc_mnt +...
           2.*mat_loc_v_1 + (4-pt).*mat_loc_v_2;
   else
       warning ('geopdes:jacdet_zero_at_quad_node',...
           'op_bordo_plate_grad: singular map in element number %d', iel)
   end
end

end