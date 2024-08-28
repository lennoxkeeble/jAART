module CriticalCurve

nanmin(array::AbstractArray) = minimum(x for x in array if !isnan(x))
nanargmin(array::AbstractArray) = argmin(isnan(x) ? Inf : x for x in array)
nanabsargmin(array::AbstractArray) = argmin(isnan(x) ? Inf : abs(x) for x in array)

nanmax(array::AbstractArray) = maximum(x for x in array if !isnan(x))
nanargmax(array::AbstractArray) = argmax(isnan(x) ? -Inf : x for x in array)
nanabsargmax(array::AbstractArray) = argmax(isnan(x) ? Inf : abs(x) for x in array)

lambda(r::Number, a::Number)::Number = (r^2 * (r - 3.0) + a^2 * (r + 1.0)) / a / (1.0 - r)
    
eta(r::Number, a::Number)::Number = r^3 / a^2 * (-r + 4.0 * (a^2 + r * (r - 2.0)) / (r - 1.0)^2)

r_critc_plus(a::Number)::Number = 2.0 * (1.0 + cos(2.0 * acos(a) / 3.0))

r_critc_minus(a::Number)::Number = 2.0 * (1.0 + cos(2.0 * acos(-a) / 3.0))
    
alpha(λ::Number, i::Number)::Number = -λ / sin(i)

betaSquared(λ::Number, η::Number, i::Number, a::Number)::Number = η + a^2 * cos(i)^2 - λ^2 * tan(i)^-2

@views function compute_critical_curve(spin::Number, inc_deg::Number, nPoints::Int64)
    α = zeros(nPoints); β = zeros(nPoints);
    inc_rad = deg2rad(inc_deg)

    r_plus = r_critc_plus(spin)
    r_minus = r_critc_minus(spin)
    r = range(start = r_minus, stop = r_plus, length = nPoints) |> collect
    @inbounds for i in eachindex(r)
        λ = lambda(r[i], spin)
        η = eta(r[i], spin)
        α[i] = alpha(λ, inc_rad)
        β[i] = betaSquared(λ, η, inc_rad, spin)

        β[i] > 0.0 ? β[i] = sqrt(β[i]) : (β[i] = NaN; α[i] = NaN)
    end
    # return filter!(!isnan, α), filter!(!isnan, β)
    return α, β
end

function rotate_x_y_points(x::Number, y::Number, angle::Number)
    return x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle)
end

function rotate_x_y_arrays!(x::AbstractArray, y::AbstractArray, x_rot::AbstractArray, y_rot::AbstractArray, angle::Number)
    for i in eachindex(x)
        x_rot[i], y_rot[i] = rotate_x_y_points(x[i], y[i], angle)
    end
end

# compute diameter at angle varphi
function compute_diameter(alpha::AbstractArray, beta::AbstractArray, alpha_rot::AbstractArray, beta_rot::AbstractArray, varphi::Number)
    # Convert varphi to radians
    varphi_rad = deg2rad(varphi)

    rotate_x_y_arrays!(alpha, beta, alpha_rot, beta_rot, -varphi_rad)
    
    # Calculate diameter (distance between positive and negative x-intercepts)
    return nanmax(alpha_rot) - nanmin(alpha_rot)
end

end