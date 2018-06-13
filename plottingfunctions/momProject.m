function [parallel_proj, perpendicular_proj] = ...
    momProject(momX, momY, momZ, mass, p, p_on)
    
    CMCX = (momX(:, p_on(1))+momX(:, p_on(2)))/sum(mass(p_on));
    CMCY = (momY(:, p_on(1))+momY(:, p_on(2)))/sum(mass(p_on));
    CMCZ = (momZ(:, p_on(1))+momZ(:, p_on(2)))/sum(mass(p_on));

    momX(:, p) = momX(:, p) - CMCX*mass(p);
    momX(:, p_on(1)) = momX(:, p_on(1)) - CMCX*mass(p_on(1));
    momX(:, p_on(2)) = momX(:, p_on(2)) - CMCX*mass(p_on(2));

    momY(:, p) = momY(:, p) - CMCX*mass(p);
    momY(:, p_on(1)) = momY(:, p_on(1)) - CMCX*mass(p_on(1));
    momY(:, p_on(2)) = momY(:, p_on(2)) - CMCX*mass(p_on(2));

    momZ(:, p) = momZ(:, p) - CMCX*mass(p);
    momZ(:, p_on(1)) = momZ(:, p_on(1)) - CMCX*mass(p_on(1));
    momZ(:, p_on(2)) = momZ(:, p_on(2)) - CMCX*mass(p_on(2));

    V = [(momX(:, p_on(2))-momX(:, p_on(1))), ...
         (momY(:, p_on(2))-momY(:, p_on(1))), ...
         (momZ(:, p_on(2))-momZ(:, p_on(1)))];
    V = V./apply_to_rows(@norm, V);

    parallel_proj  = dot(V, [momX(:, p), momY(:, p), momZ(:, p)], 2);
    perpendicular_proj = apply_to_rows(@norm, [momX(:, p), momY(:, p), momZ(:, p)] - V.*parallel_proj); 

end