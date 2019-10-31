
using Distributions, Distributed, DelimitedFiles;

include("function.jl")

Ne = 50
r = ARGS[1]
coeffSel = ARGS[2]
u = ARGS[3]

r = parse(Float64, r)
coeffSel = parse(Float64, coeffSel)
u = parse(Float64, u)

## pdf function
# function to calculate pdf 
# per column in A matrix
function CalcPdf(; Ne, Pmut = [0.25 0.25 0.25 0.25], Xmatrix, Coeff)
    # Pmut: a 1x4 matrix used to calculate one colume in the A matrix
    k = Coeff .* prod(Pmut .^ Xmatrix, dims = 2)
    return k
end  

function Pcalc_A(; Ne, Coeff, Pmut)
    # function to calculate A Matrix
    n = 2*Ne
    nrow = Int((2*Ne+1)*(2*Ne+2)*(2*Ne+3)/6)
    X = Xcalc(Ne)
    A = Array{Float64,2}(undef, nrow, nrow)


    for i = 1:nrow
        Pmutt = Pmut[:,i]
        A[:,i] = CalcPdf(Ne = Ne, Pmut = Pmutt', Xmatrix = X, Coeff = Coeff)
    end
    return A 
end

function run_generations(Ne,nGenerations,p_init, A, interval)
    count = interval
    p = p_init
    for i = 1:nGenerations
        if count == interval 
        writedlm(output, p) 
            count = 0
        end
            count += 1
        
        # calculate updated p vector for next generation
        p = A * p
    end
    return p
end     

### Calculate the coefficient of multinomial distribution
p = initialP = initialPcal(Ne = Ne) # initial p vector
Coeff = initialP ./ 0.25^(2*Ne); # multinomial coefficients for combinations in X

# a nrowx4  matrix and each row stands for updated haplotype frequencies coresponding to each row of X matrix 
PmutG = Pmut_calc(Precomb = Precomb_calc(Ne = Ne, r = r)[1], u = u)
PselG = Psel_calc(PmutG, coeffSel)
# doesn't change from generation to generation

rsqr = Precomb_calc(Ne = Ne, r = r)[2];

## Generation and Graph 
nGenerations = 3500 # number of generations of random mating
interval = 50; # output frequency

# A matrix
@time A = Pcalc_A(Ne = Ne, Coeff = Coeff, Pmut = PselG');

filename = "selection_"*"Ne"*string(Ne)*"_r"*string(r)*"_coeffSel"*string(coeffSel)*".txt"

output = open(filename, "w")

@time run_generations(Ne, nGenerations, p, A, interval)

close(output)
