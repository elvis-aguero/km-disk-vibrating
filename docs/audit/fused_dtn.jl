# Prototype: allocation-free far_part! (same math as src/dtn/dtn_generator.jl), checked
# for bitwise-ish equality against the shipped implementation.
using FaradayDisk, Printf, Logging
global_logger(Logging.NullLogger())

nr0, D0 = 50, 20.0
ref = generate_dtn(nr0, D0)        # shipped implementation
t_ref = @elapsed generate_dtn(nr0, D0)
a_ref = @allocated generate_dtn(nr0, D0)

@eval FaradayDisk function far_part!(Line::Vector{Float64}, rn_k::Float64, i_far::Vector{Float64}, Kern_far::Vector{Float64},
                    cos_l::Vector{Float64}, sin_l::Vector{Float64}, nr::Int, refp::Int)
    @inbounds for l in eachindex(cos_l)
        c = cos_l[l]; s = sin_l[l]
        for a in eachindex(i_far)
            i = i_far[a]; K = Kern_far[a]
            re = rn_k + i * c / refp
            im = i * s / refp
            radn = abs(sqrt(re^2 + im^2))
            idx = floor(radn)
            w1 = clamp(radn - idx, 0.0, 1.0); w2 = w1^2; w3 = w1^3
            ix = Int(idx)
            if idx < 0.5
                Line[ix + 1] += -(3 * w3 / 4 - 7 * w2 / 4 + 1) * K
                Line[ix + 2] += -(-w3 + 2 * w2) * K
                Line[ix + 3] += -(w3 - w2) / 4 * K
            end
            if idx >= 0.5 && idx < nr
                Line[ix] += -(-w3 / 6 + w2 / 2 - w1 / 3) * K
                Line[ix + 1] += -(w3 / 2 - w2 - w1 / 2 + 1) * K
            end
            if idx < nr - 1 && idx >= 0.5
                Line[ix + 2] += -(-w3 / 2 + w2 / 2 + w1) * K
            end
            if idx < nr - 2 && idx >= 0.5
                Line[ix + 3] += -(w3 - w1) / 6 * K
            end
        end
    end
    return Line
end

new = generate_dtn(nr0, D0)
t_new = @elapsed generate_dtn(nr0, D0)
a_new = @allocated generate_dtn(nr0, D0)
@printf("nr=50: rel diff %.2e | shipped %.2fs %.2f GB | fused %.3fs %.4f GB | speedup %.1fx\n",
        maximum(abs, new - ref) / maximum(abs, ref), t_ref, a_ref / 1e9, t_new, a_new / 1e9, t_ref / t_new)
for nr in (100, 200)
    t = @elapsed generate_dtn(nr, 0.4nr)
    @printf("fused nr=%d: %.2fs\n", nr, t)
end
