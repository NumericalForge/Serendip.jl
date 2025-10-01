


function set_joint_rotation(J::AbstractMatrix{Float64}, R::AbstractMatrix{Float64})
    if size(J,2) == 2
        # 3D interface: J ∈ ℝ^{3×2}, columns are tangents v2, v3
        @inbounds begin
            v2x, v2y, v2z = J[:,1]
            v3x, v3y, v3z = J[:,2]

            # v1 = v2 × v3
            r1x = v2y*v3z - v2z*v3y
            r1y = v2z*v3x - v2x*v3z
            r1z = v2x*v3y - v2y*v3x

            # re-orthogonalize v2 := v3 × v1
            r2x = v3y*r1z - v3z*r1y
            r2y = v3z*r1x - v3x*r1z
            r2z = v3x*r1y - v3y*r1x

            # v3 stays as given
            r3x, r3y, r3z = v3x, v3y, v3z

            # normalize rows
            n1 = sqrt(r1x*r1x + r1y*r1y + r1z*r1z)
            n2 = sqrt(r2x*r2x + r2y*r2y + r2z*r2z)
            n3 = sqrt(r3x*r3x + r3y*r3y + r3z*r3z)

            # write into R (3×3)
            R[:, 1] .= (r1x/n1, r1y/n1, r1z/n1)
            R[:, 2] .= (r2x/n2, r2y/n2, r2z/n2)
            R[:, 3] .= (r3x/n3, r3y/n3, r3z/n3)
        end
    else
        # 2D interface: J ∈ ℝ^{2×1}, column is tangent v2
        @inbounds begin
            v2x, v2y = J[1,1], J[2,1]

            # v1 ⟂ v2 (counter-clockwise)
            r1x, r1y = v2y, -v2x
            r2x, r2y = v2x,  v2y

            # normalize rows
            n1 = sqrt(r1x*r1x + r1y*r1y)
            n2 = sqrt(r2x*r2x + r2y*r2y)

            R[:, 1] .= (r1x/n1, r1y/n1)
            R[:, 2] .= (r2x/n2, r2y/n2)
        end
    end
    return R
end
