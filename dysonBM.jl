using Random
using LinearAlgebra
using DataFrames, CSV

"""
    genGUEProcess(m, n, rng)
    Generates `n` step GUE process with matrix size `m`.
    `rng` : random number generator
"""
function genGUEProcess(m, n, rng)
    matrixLs = Matrix{ComplexF16}[]
    
    #initialize
    for t in 1:n+1
        push!(matrixLs, zeros(ComplexF16, m, m))
    end

    #sampling
    for i in 1:m
        for j in i:m
            if i == j
                motion = brownianMotion(n, rng)
                for t in 1:n+1
                    matrixLs[t][i,j] = motion[t]
                end
            else
                motionR = brownianMotion(n, rng)
                motionI = brownianMotion(n, rng)
                for t in 1:n+1
                    matrixLs[t][i,j] = (motionR[t] + motionI[t] * im)/sqrt(2)
                    matrixLs[t][j,i] = (motionR[t] - motionI[t] * im)/sqrt(2)
                end
            end
        end
    end

    return matrixLs
end

"""
    brownianMotion(n,rng)
    Generate `n` step standard Brownian motion as scaling of random walk
    `rng` : random number generator
"""
function brownianMotion(n, rng)
    motion = zeros(Float16, n+1)
    
    for step in 1:n
        motion[step+1] = motion[step] + rand(rng, [-1,1])/sqrt(n)
    end

    return motion
end

"""
    eigvalsMatProcess(s)
    Returns the eigenvalue process of a matrix process `s`
    `s`::Vector{Matrix{ComplexF16}}
"""
function eigvalsMatProcess(s)
    n = size(s,1)
    eigenLs = Vector{Float16}[]

    for t in 1:n
        push!(eigenLs,eigvals(s[t]))
    end

    return eigenLs
end

"""
"""
function exportData(A, fileName)
    df = DataFrame(A, :auto)
    CSV.write(fileName,df)
end

rng = MersenneTwister(4649)
m = 5
n = 500

matrixLs = genGUEProcess(m, n, rng)
eigenLs = eigvalsMatProcess(matrixLs)
exportData(eigenLs,"DysonTest.csv")