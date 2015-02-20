%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% CALCOLA LE PARTI DI BORDO CHE NON VANNO MAI A ZERO                      %
% PER IL PROBLEMA ANTI-PLANE.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = op_Bordo_Fracture (spu, spv, msh)

mat = spalloc(spu.ndof, spu.ndof, 1);

spv.connectivity = spu.connectivity;

%keyboard
for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
        %mat_loc = zeros ( spu.nsh(iel), spu.nsh(iel) );
        mat_loc1 = zeros ( spu.nsh(iel), spu.nsh(iel) );
        mat_loc2 = zeros ( spu.nsh(iel), spu.nsh(iel) );
        
        for idof = 1:spu.nsh(iel)
            %ispv = reshape (spu.shape_functions(:,idof,iel),[], msh.nqn); % shape function v
            isnv = spv.deriv_normal(:,:,idof,iel);                        % shape deriv normal
            for jdof = 1:spu.nsh(iel)
                jsnnu =  spv.deriv2_normal(:,:,jdof,iel);      % shape deriv2 normal
                %jsnnnu = spv.deriv3_normal(:,:,jdof,iel);      % shape deriv3 normal
                % The cycle on the quadrature points is vectorized
                %for inode = 1:msh.nqn
                % TRACTION PART
                %mat_loc1(idof, jdof) = mat_loc1(idof, jdof) + ...
                %    sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*sum ( jsnnnu(3,:) .* ispv,1).');
                % COUPLE PART
                mat_loc2(idof, jdof) = mat_loc2(idof, jdof) + ...
                    sum (msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*sum ( jsnnu(3,:).*isnv(1,:),1).'); 
                
                %mat_loc = mat_loc + mat_loc1 + mat_loc2;
                %end
            end
        end
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
            mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) +  mat_loc1 + mat_loc2;
    else
        warning ('geopdes:jacdet_zero_at_quad_node', 'op_Bordo_Couple: singular map in element number %d', iel)
    end
    
end

end