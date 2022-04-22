using .problemStructsAndFunctions: JacobianPiece, ProblemConditions

function retrievePreviousConditions(current_conditions::ProblemConditions, method::String)::ProblemConditions
    return deepcopy(current_conditions)
end
