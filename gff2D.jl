using Random
using DataFrames, CSV

"""
    genGFF2d(m, rng)
    Generate a 2D GFF on square with `m` cutoff of momenta
"""
function genGFF2d(m, rng)
    s = zeros(Float16, m+1, m+1)

    for k in 0:m
        for l in 0:m
            if k == 0 && l == 0
                a = 0.0
            else
                a = randn(rng, Float16)/sqrt(k^2 + l^2)
            end

            for i in 1:m+1
                for j in 1:m+1
                    s[i,j] += a * sin(π * (i-1) * k / m) * sin(π * (j-1) * l /m)
                end
            end
        end
    end

    return s
end

"""
"""
function exportData(A, fileName)
    df = DataFrame(A, :auto)
    CSV.write(fileName,df)
end

m = 200
rng = MersenneTwister(4649)

s = genGFF2d(m, rng)
display(s)

exportData(s,"GFFtest.csv")