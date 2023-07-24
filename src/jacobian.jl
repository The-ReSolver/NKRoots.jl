export jacobianFD,
       jacobianAD,
       update!

using DualNumbers
 
struct JacobianAD{X, FT, TUP}
    F::FT  # operator
    x::X   # point at which the derivative is taken
  tmp::TUP # temporaries with dual
end

# outer constructor
jacobianAD(F, x) = 
    JacobianAD(F, x, (similar(x, Dual{eltype(x)}), similar(x, Dual{eltype(x)})))

# does nothing
update!(j::JacobianAD) = nothing

# matrix vector product function required by GMRES
function Base.:*(op::JacobianAD{X}, y::X) where {X}
    # create output
    out = similar(y)

    # aliases
     x_dual = op.tmp[1]
    Fx_dual = op.tmp[2]

    # perturbed input
    x_dual .= op.x .+ DualNumbers.ε .* y
    
    # calculate function with perturbed input
    op.F(Fx_dual, x_dual)

    # exact gradient 
    out .= dualpart.(Fx_dual)

    return out
end

# Object to approximate the action of the jacobian of a function
# using a finite difference quotient.
# Assumes F obeys :
# F(out, x)
# and overwrites `out` with F(x), returning out
struct JacobianFD{X, FT}
      F::FT          # operator
      x::X           # point at which the derivative is taken
     Fx::X           # value of the function at `x`
    tmp::Tuple{X, X} # temporaries
      ϵ::Float64     # finite difference step
end

# outer constructor
jacobianFD(F, x, ϵ::Float64=1e-6) =
    JacobianFD(F, x, F(similar(x), x), (similar(x), similar(x)), ϵ)

# update Fx assuming x has been modified outside
function update!(j::JacobianFD{X}) where {X}
    j.F(j.Fx, j.x)
    return nothing
end

# matrix vector product function required by GMRES
function Base.:*(op::JacobianFD{X}, y::X) where {X}
    # create output
    out = similar(y)

    # aliases
     x_plus_y = op.tmp[1]
    Fx_plus_y = op.tmp[2]
            x = op.x
           Fx = op.Fx
            ϵ = op.ϵ
            δ = ϵ.*norm(y)
    
    # if perturbation is zero set output to zero
    if norm(y) == 0
        out .= 0
    else
        # perturbed input
        x_plus_y .= x .+ δ .* y
        
        # calculate function with perturbed input
        op.F(Fx_plus_y, x_plus_y)

        # calculate finite difference approximation of gradient 
        out .= (Fx_plus_y .- Fx)./δ
    end

    return out
end

