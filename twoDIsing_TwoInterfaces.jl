using Random
using DataFrames, CSV

"""
    ising2d_twoInterfaces!(s, β, niters, rng)
    Update configurations of 2d Ising, boundary conditions change at the four corners
"""
function ising2d_twoInterfaces!(s, β, niters, rng)
    m, n = size(s)
    min_h = -4
    max_h = 4
    prob = [1/(1+exp(-2*β*h)) for h in min_h:max_h]
    for iter in 1:niters
        for j in 1:n
            for i in 1:m
                if (i==1)
                    NN = 1
                else
                    NN = s[i-1,j]
                end

                if (i==m)
                    SS = 1
                else
                    SS = s[i+1,j]
                end

                if (j==1)
                    WW = -1
                else
                    WW = s[i,j-1]
                end
                
                if (j==n)
                    EE = -1
                else
                    EE = s[i,j+1]
                end
                
                h = NN + SS + WW + EE
                s[i,j] = ifelse(rand(rng) < prob[h-min_h+1], +1, -1)
            end
        end
    end
end

"""
"""
function exportData(A, fileName)
    df = DataFrame(A, :auto)
    CSV.write(fileName,df)
end

const β_crit = log(1+sqrt(2))/2
const m = 200
rng = MersenneTwister(5912)
s = rand(rng, Int8[-1,1], m, m)

ising2d_twoInterfaces!(s, β_crit, 4000, rng)

exportData(s,"IsingTest.csv")