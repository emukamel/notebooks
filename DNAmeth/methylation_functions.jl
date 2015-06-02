using PyPlot
import Color: convert, HSL, RGB
using Optim

function get_colors(N)
    h = 0
    col = zeros(3,N)
    for i = 1:N
        color = convert(RGB, HSL((h%1)*360, 0.5, 0.5));
        col[:,i] = [color.r color.g color.b]
        h += 0.618033988749895 # silver ratio
    end
    return col
end

function plot_meth(X,c)
    nc = int(maximum(c)) # number of clusters
    if nc > 1
        col = get_colors(nc)
    else
        col = zeros(3,1) # black
    end
    N,W = size(X)
    ind = sortperm(c)
    y = 1
    for n = ind
        plot([0.5,W+0.5],[y,y],"-",color=col[:,c[n]])
        for w = 1:W
        	if isnan(X[n,w])
        		continue
	        elseif X[n,w] == 1.0
                plot(w,y,"o",color=col[:,c[n]])
            else
                plot(w,y,"o",color=col[:,c[n]],mfc="w")
            end
        end
        y += 1
    end
    xlim([0,1+W])
    ylim([0,1+N])
    xticks([])
    yticks([])
    if nc > 1
        title("All reads, with assigned clusters")
    else
        title("All reads")
    end
end

function plot_meth_avg(X,c)
    nc = int(maximum(c)) # number of clusters
    col = get_colors(nc)
    N,W = size(X)
    
    for cc = 1:nc
        plot([0.5,W+0.5],[cc,cc],"-",color=col[:,cc])
        Xavg = zeros(1,W)
        for w = 1:W
            xx = (1+X[c .== cc,w])/2
            xavg = mean(xx[!isnan(xx)])
            color = xavg*col[:,cc] + (1-xavg)*[1,1,1];
            plot(w,cc,"o",color=color)
        end
    end
    xlim([0,1+W])
    ylim([0,1+nc])
    xticks([])
    yticks([])
    title("average methylation pattern")
end

function blackwell_macqueen(;N=30,W=10)
    # Blackwell-MacQueen urn scheme to generate reads
    # -----------------------------
    # N = number of reads
    # W = methylation sites per read

    X = bool(zeros(N,W))
    P = zeros(1,W)
    G = {"a"=>ones(1,W),"b"=>ones(1,W)} # alpha and beta for beta prior
    alpha = 1 # dispersion parameter for Dirichlet Process

    # initialize probabilities for first cluster
    for w = 1:W
        a,b = G["a"][w],G["b"][w]
        P[1,w] = rand(Beta(a,b))
    end
    c = zeros(N)
    c[1] = 1

    for n = 2:N
        # sample from categorical distribution with probability vector p_
        #  -- A new cluster is formed with probability alpha/(alpha+n-1)
        #  -- Each old datapoint has probability 1/(alpha+n-1) of being reused
        p_ = [ones(n-1)/(n-1+alpha), alpha/(n-1+alpha)]
        j = rand(Categorical(p_))
        if j <= (n-1) 
            c[n] = c[j] # Re-use cluster c[j]
        else
            # Make a new cluster
            P = [P;zeros(1,W)]
            for w = 1:W
                a,b = G["a"][w],G["b"][w]
                P[end,w] = rand(Beta(a,b))
            end
            c[n] = size(P,1) # new cluster index
        end

        for w = 1:W
            X[n,w] = rand(Bernoulli(P[c[n],w]))
        end
    end
    
    return X,P,c
end

function beta_binomial_likelihood(Xc,G)
    # Computes the likelihood for a subset of the data, Xc
    n = size(Xc,1)
    x_sum = sum(Xc,1)
    return sum(lgamma(n + 1) + 
               lgamma(x_sum + G["a"]) + 
               lgamma(n - x_sum + G["b"]) + 
               lgamma(G["a"] + G["b"]) - 
               lgamma(x_sum + 1) - 
               lgamma(n - x_sum + 1) - 
               lgamma(n + G["a"] + G["b"]) - 
               lgamma(G["a"]) - 
               lgamma(G["b"]))
end

function generate_reads(;N=100,W=21,L=100,Cprob=[0.3, 0.7],alpha=1,beta=1,targ_dist=NaN)
	# Generate synthetic data of overlapping reads
    # -----------------------------
    # N = number of reads
    # W = methylation sites per read
    # L = number of methylation-eligible sites
    # Cprob = array listing the proportions of cell types
    # alpha/beta = parameters for beta distribution (prior)
    X = fill(NaN,(N,L))
    c = zeros(N)
    nC = length(Cprob) # number of cell types
    cumC = cumsum(Cprob) - Cprob[1]

    if isnan(targ_dist)
        # generate random methylation probabilities
        P = rand(Beta(alpha,beta),(nC,L))
    else
        # generate probabilities that have specified distance
        assert(nC == 2) 
        P = generate_patterns_JSdist(L,targ_dist)
    end

    if floor(W/2) == W/2
    	W += 1 # make sure W is odd
    end

    for i = 1:N
    	m = rand(1:L) # middle point of read
    	c[i] = find((i/N) .> cumC)[end] # current cluster
    	s1 = maximum([1,m-floor(W/2)])
    	s2 = minimum([L,m+floor(W/2)])

    	for w = s1:s2
    		if rand(Bernoulli(P[c[i],w])) == 1
    			X[i,w] = 1.0
    		else
    			X[i,w] = -1.0
    		end
    	end
    end

    return X,P,c
end



##
## Functions to initialize methylation probabilities
##

function d_kl_bern(p,q)
    ## Kulback-Leibler divergence for two Bernoulli variables with probabilities p and q
    return p*log2(p./q)+(1-p)*log2((1-p)/(1-q))
end

function d_js_bern(p1,p2)
    ## Jenson-Shannon divergence for two Bernoulli variables with probabilities p1 and p2
    m = 0.5*(p1+p2)
    return 0.5*d_kl_bern(p1,m) + 0.5*d_kl_bern(p2,m)
end

function d_js_multivar_bern(x::Vector,y::Vector)
    ## Jenson-Shannon divergence between two multivariate Bernoulli variables
    dist = 0
    for i=1:length(x)/2
        dist += d_js_bern(x[i],y[i])
    end
    return dist
end

function impose_bounds!(x::Vector)
    x[x.>1] = 1
    x[x.<=0] = eps()
end

function generate_patterns_JSdist(n,targ_dist;max_iter=20)
    ## n = number of methylation sites
    ## targ_dist = target Jensen-Shannon divergence
    function f(x::Vector)
        impose_bounds!(x)
        return abs(d_js_multivar_bern(x[1:end/2],x[1+end/2:end]))
    end
    df = Optim.autodiff(x->abs(f(x)-targ_dist), Float64, n*2)
    
    converged = false
    iter = 0
    y = NaN
    while !converged && iter < max_iter
        try
            y = fminbox(df,rand(n*2),zeros(n*2),ones(n*2))
            converged = true
        catch
            println("didnt converge, trying again")
            iter +=1
        end
    end
    if converged
        p1,p2 = y.minimum[1:end/2],y.minimum[1+end/2:end]
        return [p1'; p2']
    else
        error("Unable to find solution")
    end
end