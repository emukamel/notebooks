function run_fit(A,k,loss,rx,ry,initX,initY)
    ## Fits a low-rank model given:
    ##   loss()      -- the loss function
    ##   rx(),ry()   -- regularization functions
    ##   initX,initY -- initial guesses for X,Y
    losses = fill(loss,size(A,2))
    obs = find_observations(A);
    glrm = GLRM(A,obs,losses,rx,ry,k)

    if initX != nothing
        glrm.X = initX
    end
    if initY != nothing
        glrm.Y = initY
    end
    X,Y,ch = fit!(glrm,verbose=false)
    return X,Y,ch
end

function run_batch(A,k,loss_function,rx,ry,initX,initY,N)
    ## Fits a specified low-rank model N times, returns best one
    best = Inf
    X,Y,ch = 0,0,0
    for i = 1:N
        Xest,Yest,ch_ = run_fit(A,k,loss_function,rx,ry,initX,initY)
        if ch_.objective[end] < best
            best = ch_.objective[end]
            X,Y,ch = Xest,Yest,ch_
        end
    end
    println("Best Objective: ",ch.objective[end])
    return X,Y,ch
end

function find_observations(A)
    m,n = size(A)
    obs = (Int64,Int64)[]
    for i = 1:m
        for j = 1:n
            if ~isnan(A[i,j])
                push!(obs,(i,j))
            end
        end
    end
    return obs
end