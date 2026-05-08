# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Constraint{T}
    kind    ::Symbol
    selector::Any
    expr    ::Expr
    target  ::Vector{T}
    terms   ::Vector{Pair{Symbol,Float64}}
    rhs     ::Float64
end


function _constraint_target_nodes(constraint::Constraint)
    if constraint.kind == :node
        return unique(node -> node.id, constraint.target)
    end
    return unique(node -> node.id, Node[node for facet in constraint.target for node in facet.nodes])
end


function mount_constraint_matrices(constraints::Vector{Constraint}, nu::Int, ndofs::Int)
    # First assemble the full constraint system A*u = b using the global dof
    # numbering. In the current static LM implementation, every constraint row
    # must depend only on unknown dofs; constraints touching prescribed dofs
    # are rejected to keep the reduced solve and residual bookkeeping simple.
    R, C, V = Int[], Int[], Float64[]
    b_full = Float64[]

    for constraint in constraints
        # Each constraint expression generates one algebraic row for every
        # selected node.
        nodes = _constraint_target_nodes(constraint)
        for node in nodes
            row = length(b_full) + 1
            push!(b_full, constraint.rhs)

            for (key, coeff) in constraint.terms
                dof = get_dof(node, key)
                dof === nothing && error("mount_constraint_matrices: Node $(node.id) has no dof `$(key)` required by constraint $(constraint.expr)")
                push!(R, row)
                push!(C, dof.eq_id)
                push!(V, coeff)
            end
        end
    end

    A = sparse(R, C, V, length(b_full), ndofs)

    # Reject rows that touch prescribed dofs. Such mixed constraints are
    # mathematically valid, but they are intentionally out of scope for this
    # first implementation.
    if nu < ndofs
        for row in 1:size(A, 1)
            nnz(A[row, nu+1:ndofs]) == 0 && continue
            error("mount_constraint_matrices: Constraints involving prescribed dofs are not supported in this version.")
        end
    end

    # With the restriction above, the full rows are already the active LM rows.
    A1 = A[:, 1:nu]
    b = b_full

    return A1, b
end
