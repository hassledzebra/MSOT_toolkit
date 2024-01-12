% Joins two structures together to make another structure. structB
% dominates structA, so if structA and structB have the same fields,
% structB's value will be the one in fusedStructs.
function fusedStructs = fuseStructures(structA,structB)
   
    fusedStructs = structA;
    
    Bfields = fieldnames(structB);
    NBfields = numel(Bfields);
    
    for k = 1:NBfields
        fusedStructs.(Bfields{k}) = structB.(Bfields{k});
    end
    
end