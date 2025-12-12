# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Tensor definitions using Mandel notation

const I2 = SVector(1., 1., 1., 0., 0., 0.)
const I4 = SMatrix{6,6}(I)

const Psd = @SArray [
    2/3. -1/3. -1/3. 0. 0. 0.
   -1/3.  2/3. -1/3. 0. 0. 0.
   -1/3. -1/3.  2/3. 0. 0. 0.
      0.    0.    0. 1. 0. 0.
      0.    0.    0. 0. 1. 0.
      0.    0.    0. 0. 0. 1. ]

dev(T::Vec6) = Psd*T # deviatoric tensor

# Tensor invariants
LinearAlgebra.tr(σ::Vec6) = σ[1]+σ[2]+σ[3]

# Deviatoric tensor invariants (Mandel notation)
function J2(σ::Vec6)
    t11, t22, t33, t23, t13, t12 = σ[1], σ[2], σ[3], σ[4]/SR2, σ[5]/SR2, σ[6]/SR2
    return 1/6*( (t11-t22)^2 + (t22-t33)^2 + (t33-t11)^2 ) + t23^2 + t13^2 + t12^2
end


# Third invariant of the deviatoric tensor
function J3(σ::Vec6)
    t11, t22, t33 = σ[1], σ[2], σ[3]
    s23, s13, s12 = σ[4]/SR2, σ[5]/SR2, σ[6]/SR2
    p = t11+t22+t33
    s11 = t11 - 1/3*p
    s22 = t22 - 1/3*p
    s33 = t33 - 1/3*p
    return s11*s22*s33 + 2*s12*s23*s13 - s11*s23^2 - s22*s13^2 - s33*s12^2
end


# """
# This function is not precise enouugh...
# """
# function eigvals2(T::Vec6)
#     @assert length(T) == 6

#     t11, t22, t33, t12, t23, t13 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2

#     i1 = t11 + t22 + t33
#     i2 = t11*t22 + t22*t33 + t11*t33 - t12*t12 - t23*t23 - t13*t13

#     i1==0.0 && i2==0.0 && return zeros(3)

#     i3  = t11*(t22*t33 - t23*t23) - t12*(t12*t33 - t23*t13) + t13*(t12*t23 - t22*t13)
#     val = (2*i1^3 - 9*i1*i2 + 27*i3 )/( 2*(i1^2 - 3*i2)^(3/2) )
#     val = clamp(val, -1.0, 1.0) # to avoid 1.000000000000001

#     θ = 1/3*acos( val )

#     r = 2/3*√(i1^2-3*i2)

#     s1 = i1/3 + r*cos(θ)
#     s2 = i1/3 + r*cos(θ - 2*π/3)
#     s3 = i1/3 + r*cos(θ - 4*π/3)

#     # sorting
#     if s1<s2; s1,s2 = s2,s1 end
#     if s2<s3; s2,s3 = s3,s2 end
#     if s1<s2; s1,s2 = s2,s1 end

#     return Vec3(s1, s2, s3)
# end



# function eigvals(i1::Float64, j2::Float64, θ::Float64)
#     # 1. Hydrostatic component (center of the Mohr circle)
#     p = i1 * 0.3333333333333333  # 1/3

#     # 2. Deviatoric Radius
#     #    If j2 is effectively zero, the stress is purely hydrostatic
#     if j2 <= 1e-16
#         return Vec3(p, p, p)
#     end
    
#     # r = 2/sqrt(3) * sqrt(J2)
#     r = 1.1547005383792515 * sqrt(j2)

#     # 3. Deviatoric Eigenvalues (based on Lode angle)
#     #    Using the standard trigonometric solution
#     #    Note: 2π/3 ≈ 2.0943951023931953
#     s1 = r * cos(θ)
#     s2 = r * cos(θ - 2.0943951023931953)
#     s3 = r * cos(θ + 2.0943951023931953)

#     # 4. Total Eigenvalues
#     λ1 = p + s1
#     λ2 = p + s2
#     λ3 = p + s3

#     # 5. Sort Descending (λ1 >= λ2 >= λ3)
#     #    Manual sort is faster than allocating a vector for sort()
#     if λ1 < λ2; λ1, λ2 = λ2, λ1; end
#     if λ2 < λ3; λ2, λ3 = λ3, λ2; end
#     if λ1 < λ2; λ1, λ2 = λ2, λ1; end

#     return Vec3(λ1, λ2, λ3)
# end



# """
# Computes the eigenvalues of a second order tensor written in Mandel notation.
# The eigenvalues are sorted from highest to lowest
# """
# function eigvals(T::Vec6; sort=true)
#     t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2

#     # full notation
#     F = @SArray[ t11  t12  t13
#                  t12  t22  t23
#                  t13  t23  t33 ]

#     L, _ = eigen(F, permute=false, scale=false)

#     # put biggest eigenvalue first
#     sort && return Base.sort(L, rev=true)
#     return L
# end


"""
    eigenvalues(T::Vec6)

Compute the eigenvalues of a symmetric tensor T (in Voigt notation) using 
the analytical method. 
Internally uses deviatoric invariants (J2, J3) for maximum numerical stability.
"""
function eigvals(T::Vec6)
    # 1. Extract components (handle Mandel notation if necessary)
    #    Assuming T = [11, 22, 33, 23, 13, 12] or similar.
    #    Adjust the sqrt(2) factors if your Vec6 is Mandel.
    t11, t22, t33 = T[1], T[2], T[3]
    t23, t13, t12 = T[4]/SR2, T[5]/SR2, T[6]/SR2 

    # 2. Compute I1 (Hydrostatic Part)
    i1 = t11 + t22 + t33
    p  = i1 * 0.3333333333333333 # Mean normal stress

    # 3. Compute Deviatoric Stress (s = t - p*I)
    s11, s22, s33 = t11 - p, t22 - p, t33 - p
    
    # 4. Compute J2 (Deviatoric Second Invariant)
    #    Using s_ij is numerically safer than using I2
    j2 = 0.5*(s11^2 + s22^2 + s33^2) + t12^2 + t23^2 + t13^2

    # 5. Handle the Hydrostatic/Zero Case
    if j2 <= 1e-16
        # If deviatoric energy is zero, all eigenvalues = mean stress
        return Vec3(p, p, p)
    end

    # 6. Compute J3 (Deviatoric Third Invariant) for Angle
    j3 = s11*s22*s33 + 2.0*t12*t23*t13 - s11*t23^2 - s22*t13^2 - s33*t12^2

    # 7. Calculate Lode Angle θ
    #    Normalized J3 argument for acos
    arg = (1.5 * sqrt(3.0)) * j3 / (j2 * sqrt(j2))
    arg = clamp(arg, -1.0, 1.0)
    θ   = acos(arg) * 0.3333333333333333

    # 8. Solve for Eigenvalues
    r = 2.0 * sqrt(j2 * 0.3333333333333333)
    # Precompute constants
    c1 = cos(θ)
    c2 = cos(θ - 2.0943951023931953) # 2π/3
    c3 = cos(θ + 2.0943951023931953) # 4π/3

    λ1 = p + r * c1
    λ2 = p + r * c2
    λ3 = p + r * c3

    # 9. Sort Descending
    if λ1 < λ2; λ1, λ2 = λ2, λ1; end
    if λ2 < λ3; λ2, λ3 = λ3, λ2; end
    if λ1 < λ2; λ1, λ2 = λ2, λ1; end

    return Vec3(λ1, λ2, λ3)
end


"""
Computes eigenvalues and eigenvectors of a second order tensor written in Mandel notation.
The first eigenvalues corresponds to the highest.
The eigenvectors are returned columnwise and disposed in a clockwise coordinate system.
"""
function LinearAlgebra.eigen(T::Vec6)
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2

    # full notation
    F = @SArray[ t11  t12  t13
                 t12  t22  t23
                 t13  t23  t33 ]

    L, V = eigen(F, permute=false, scale=false)

    # put biggest eigenvalue first
    p = sortperm(L, rev=true)
    L = L[p]
    V = V[:, p]
    V = [ V[:,1] V[:,2] normalize(cross(V[:,1], V[:,2])) ]
    # V[:,3] = normalize(cross(V[:,1], V[:,2]))

    return L, V
end


function LinearAlgebra.inv(T::Vec6)
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2

    # full notation
    F = @SArray[ t11  t12  t13
                 t12  t22  t23
                 t13  t23  t33 ]

    G = inv(F)
    return Vec6( G[1,1], G[2,2], G[3,3], SR2*G[2,3], SR2*G[1,3], SR2*G[1,2] )
end


function dott(T::AbstractVector{Float64}, V::AbstractVector{Float64})
    @assert length(T) == 6
    @assert length(V) == 3
    t11, t22, t33, t23, t13, t12 = T[1], T[2], T[3], T[4]/SR2, T[5]/SR2, T[6]/SR2
    v1, v2, v3 = V 
    return Vec3( 
        t11*v1 + t12*v2 + t13*v3,
        t12*v1 + t22*v2 + t23*v3,
        t13*v1 + t23*v2 + t33*v3
    )

end



function rotation_tensor!(V::Array{Float64,2})
    # V : second order tensor with direction cosines (new system axes in old system coordinates)
    # R : fourth order tensor

    lx, ly, lz = V[1,:]
    mx, my, mz = V[2,:]
    nx, ny, nz = V[3,:]

    return @SArray [
            lx*lx      ly*ly      lz*lz    SR2*ly*lz    SR2*lz*lx    SR2*lx*ly
            mx*mx      my*my      mz*mz    SR2*my*mz    SR2*mz*mx    SR2*mx*my
            nx*nx      ny*ny      nz*nz    SR2*ny*nz    SR2*nz*nx    SR2*nx*ny
        SR2*mx*nx  SR2*my*ny  SR2*mz*nz  my*nz+ny*mz  mx*nz+nx*mz  mx*ny+nx*my
        SR2*nx*lx  SR2*ny*ly  SR2*nz*lz  ny*lz+ly*nz  nx*lz+lx*nz  nx*ly+lx*ny
        SR2*lx*mx  SR2*ly*my  SR2*lz*mz  ly*mz+my*lz  lx*mz+mx*lz  lx*my+mx*ly ]
end



"""
Return a dictionary with conventional stress and stress values
from stress and strain tensors defined in Mandel notation.
"""
@inline function stress_strain_dict(σ::Vec6, ε::Vec6, stress_state::Symbol)
    σvm = √(3*J2(σ))
    εv = ε[1] + ε[2] + ε[3]

    if stress_state in (:plane_stress,:plane_strain)
        s1, _, s3 = eigvals(σ)
        return OrderedDict{Symbol,Float64}(
            :σxx => σ[1],
            :σyy => σ[2],
            :σzz => σ[3],
            :σyz => σ[4]/SR2,
            :σxz => σ[5]/SR2,
            :σxy => σ[6]/SR2,
            :σvm => σvm,
            :σ1  => s1,
            :σ3  => s3,
            :εxx => ε[1],
            :εyy => ε[2],
            :εzz => ε[3],
            :εxy => ε[6]/SR2,
            :εv  => εv,
        )
    elseif stress_state==:axisymmetric
        return OrderedDict{Symbol,Float64}(
            :σrr => σ[1],
            :σyy => σ[2],
            :σtt => σ[3],
            :σry => σ[6]/SR2,
            :σvm => σvm,
            :εrr => ε[1],
            :εyy => ε[2],
            :εtt => ε[3],
            :εry => ε[6]/SR2,
            :εv  => εv,
        )
    else
        s1, s2, s3 = eigvals(σ)
        return OrderedDict{Symbol,Float64}(
            :σxx => σ[1],
            :σyy => σ[2],
            :σzz => σ[3],
            :σyz => σ[4]/SR2,
            :σxz => σ[5]/SR2,
            :σxy => σ[6]/SR2,
            :σvm => σvm,
            :σ1  => s1,
            :σ2  => s2,
            :σ3  => s3,
            :εxx => ε[1],
            :εyy => ε[2],
            :εzz => ε[3],
            :εyz => ε[4]/SR2,
            :εxz => ε[5]/SR2,
            :εxy => ε[6]/SR2,
            :εv  => εv,
        )
    end

end