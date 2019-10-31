# Construct X matrix
function Xcalc(Ne)
    n = 2*Ne
    # number of combinations in X matrix
    nrow = Int((2*Ne+1)*(2*Ne+2)*(2*Ne+3)/6)
    X = Array{Int, 2}(undef, nrow, 4)
    rowi = 1
    for i = 0:n
        for j = 0:(n-i)
            for k = 0:(n-i-j)
                l = n-i-j-k
                X[rowi, :] = [i j k l]
                rowi += 1
                end 
            end
    end
    return X
end

function Precomb_calc(; Ne, r = 0.01)
    
    # Prob changes after recombination
    n = 2*Ne
    nrow = Int((2*Ne+1)*(2*Ne+2)*(2*Ne+3)/6)
    nMinus1Inv = 1/(n-1)
    Precomb = Array{Float64, 2}(undef, nrow, 4);
    
    X = Xcalc(Ne)
    X00 = X[:,1] # count for 00 haplotype from X matrix
    X01 = X[:,2] # count for 01 haplotype from X matrix 
    X10 = X[:,3] # count for 10 haplotype from X matrix
    X11 = X[:,4] # count for 11 haplotype from X matrix
    P00 = X00./n # probability for 00 haplotype from X matrix
    P01 = X01./n # probability for 01 haplotype from X matrix
    P10 = X10./n # probability for 10 haplotype from X matrix
    P11 = X11./n # probability for 11 haplotype from X matrix
    
    P1  = P10 + P11 # vector of allele frequency of locus 1 
    P2  = P01 + P11 # vector of allele frequency of locus 2
    cov = P11 - P1.* P2 # vector of cov between the loci
    rsqr=cov.^2 ./ (P1 .* (1 .- P1) .* P2 .* (1 .- P2)) # vector of r2 
    # r2 = squared correlation between the loci
    
    ##prob. to get a haplotype given recombination
    # P'00, P'01, P'10, P'11
    Precomb[:, 1] = (1-r) .* P00 .+ r*nMinus1Inv*(P00 .* ((X00 .- 1).+ X10) .+ P01 .* (X00 .+ X10)) # 00 can come from 00 without recombination, from recombination of 00 with 00 or 10, and from recombination of 01 with 00 or 10
    Precomb[:, 2] = (1-r) .* P01 .+ r*nMinus1Inv*(P00 .* (X01 .+ X11) .+ P01 .* (X01 .- 1 .+ X11)) # 01 can come from 01 without recombination, from recombination of 00 with 01 or 11, and from recombination of 01 with 01 or 11.
    Precomb[:, 3] = (1-r) .* P10 .+ r*nMinus1Inv*(P10 .* (X10 .- 1 .+ X00) + P11 .* (X10 .+ X00)) # 10 can come from 10 without recombination, from recombination of 10 with 00 or 10, and from recombination of 11 with 10 or 00.
    Precomb[:, 4] = (1-r) .* P11 .+ r*nMinus1Inv*(P10 .* (X11 .+ X01) + P11 .*(X11 .- 1 + X01)) # 11 can come from 11 without recombination, from recombination of 10 with 01 or 11, and from recombination of 11 with 11 or 01.
    
    return Precomb, rsqr # matrix of prob. after combination, r2 for combinations in X matrix
end


# recombination and mutationn transition matrix
function Pmut_calc(; Precomb, u = 0.0025)
    # Prob of mutation
    a = (1-u)*(1-u) # no mutation at locus 1 and no mutation at locus 2
    b = (1-u)*u # mutation at one locus
    c = u*u # mutation at two loci
    
    # Mutation transition matrix
    Mu = [a b b c; # 00 can come from no mutation in 00, single mutation in 01, single mutation in 10, and 2 mutations in 11
          b a c b; # 01 can come from single mutation in 00, no mutation in 01, 2 mutations in 10, and single mutation in 11 
          b c a b; # 10 can come from single mutation in 00, 2 mutations in 01, no mutation in 10, and single mutation in 11
          c b b a] # 11 can come from 2 mutations in 00, single mutation in 01, single mutation in 10, and no mutation in 11
    
    # transition matrix after recombination and mutation for combinations in X matrix
    Pmut = Mu * Precomb' # Pmut will not change from generation to generation
end


# original prob. of each combination
function initialPcal(;Ne, p = 4)
    n = 2*Ne
    nrow = Int((2*Ne+1)*(2*Ne+2)*(2*Ne+3)/6)
    X = Xcalc(Ne)
    # probability of each combination given equal haplotype frequencies (0.25)
    # each combination is a solution of multinomial distribution
    initialPvec = pdf(Multinomial(n, p), X')
end

function Psel_calc(Pmut, coeffSel)
    Mut = Pmut'
    p00 = Mut[:, 1]
    p01 = Mut[:, 2]
    p10 = Mut[:, 3]
    p11 = Mut[:, 4]
    wbar = (1 - coeffSel) .* (p00 .+ p01) .+ p10 .+ p11
    p00_new = (1 - coeffSel) .* p00 ./ wbar
    p01_new = (1 - coeffSel) .* p01 ./ wbar
    p10_new = p10 ./ wbar
    p11_new = p11 ./ wbar
    return [p00_new p01_new p10_new p11_new]
end
