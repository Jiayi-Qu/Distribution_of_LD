
using Distributions, Distributed, DataFrames, DelimitedFiles;

include("functions.jl")

Ne = ARGS[1]
r = ARGS[2] # recombination rate
u = ARGS[3]
Ne = parse(Int, Ne) # effective population size
r= parse(Float64, r)
u = parse(Float64, u)


## pdf function 
# function to calculate pdf -- need to be updated
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

function run_generations(Ne,nGenerations,p_init, A)
    p = p_init
    count = interval
    
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
initialP = initialPcal(Ne = Ne) # initial p vector
Coeff = initialP ./ 0.25^(2*Ne); # multinomial coefficients for combinations in X

# a 4xnrow matrix and each column stands for updated haplotype frequencies coresponding to each row of X matrix 
PmutG = Pmut_calc(Precomb = Precomb_calc(Ne = Ne, r = r)[1], u = u) 
# doesn't change from generation to generation

p = initialPcal(Ne = Ne)
rsqr = Precomb_calc(Ne = Ne, r = r)[2];

## Generation and Graph 
nGenerations = 3500 # number of generations of random mating
interval = 50; # output frequency

# A matrix
A = Pcalc_A(Ne = Ne, Coeff = Coeff, Pmut = PmutG);

filename = "output_"*"Ne"*string(Ne)*"_r"*string(r)*"_u"*string(u)*".txt"
output = open(filename, "w")

@time run_generations(Ne,nGenerations,p, A)

close(output)

println("Done")
