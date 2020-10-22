"""
    interiorpoint(f, x0, c; 
                  L   = (x,λ)->(f(x) - dot(λ,c(x))),
                  ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
                  ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
                  ∇c = x->ForwardDiff.jacobian(c,x),
                  tol=1e-4, maxiter = 1000,
                  μ0 = 1.0, μfactor = 0.2,
                  xrange=[-2., 3.],
                  yrange=[-2.,6.], animate=true)

Find the minimum of function `f` subject to `c(x) >= 0` using a
primal-dual interior point method.
  
# Arguments
  
- `f` function to minimizie
- `x0` starting value. Must have c(x0) > 0
- `c` constraint function. Must return an array.
- `L   = (x,λ)->(f(x) - dot(λ,c(x)))` Lagrangian
- `∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x)` Derivative of Lagrangian wrt `x`
- `∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x)` Hessian of Lagrangian wrt `x`
- `∇c = x->ForwardDiff.jacobian(c,x)` Jacobian of constraints
- `tol` convergence tolerance
- `μ0` initial μ
- `μfactor` how much to decrease μ by
- `xrange` range of x-axis for animation
- `yrange` range of y-axis for animation
- `animate` whether to create an animation (if true requires length(x)==2)
- `verbosity` higher values result in more printed output during search. 0 for no output, any number > 0 for some.  
  
# Returns

- `(fmin, xmin, iter, info, animate)` tuple consisting of minimal function
  value, minimizer, number of iterations, and convergence info

"""
function interiorpoint(f, x0, c;
                       L   = (x,λ)->(f(x) - dot(λ,c(x))),
                       ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
                       ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
                       ∇c = x->ForwardDiff.jacobian(c,x),
                       tol=1e-4, maxiter = 1000,
                       μ0 = 1.0, μfactor = 0.2,
                       xrange=[-2., 3.],
                       yrange=[-2.,6.], animate=true, verbosity=0)
  fold = f(x0)
  xold = x0
  all(c(x0).>0) || error("interiorpoint requires a starting value that strictly satisfies all constraints")
  μ = μ0
  λ = μ./c(x0)
  xchange=Inf
  fchange=Inf
  iter = 0
  μiter = 0
  stuck=0

  animate = animate && length(x0)==2
  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    ct = contour(range(xrange[1],xrange[2], length=100), 
                range(yrange[1],yrange[2], length=100),
                 (x,y) -> log(f([x,y])))
    plot!(ct, xrange, 2.5 .- xrange) # add constraint 
    anim = Animation()
  end
  foc = [∇ₓL(xold,λ); λ.*c(xold)]
  while(iter < maxiter && ((xchange>tol) || (fchange>tol) || (stuck>0)
                           || norm(foc)>tol || μ>tol) )
    # Calculate the direction for updating x and λ
    Dc = ∇c(xold)
    cx = c(xold)
    foc = ∇ₓL(xold, λ)
    H = ∇²ₓL(xold,λ)
    Δ = [H   -Dc'; λ'*Dc  diagm(cx)] \ [-foc; μ .- cx.*λ]

    # Find a step size such that λ>=0 and c(x)>=0
    # The details here could surely be improved
    α = 1.0
    acceptedstep = false
    λold = copy(λ)
    x = copy(xold)
    while (α > 1e-10)
      x = xold + α*Δ[1:length(xold)]
      λ = λold + α*Δ[(length(xold)+1):length(Δ)]
      if (all(λ.>=0) && all(c(x).>=0))
        acceptedstep=true
        break
      end
      α *= 0.5
    end
    if !acceptedstep
      stuck = 1
      break
    end
    fnew = f(x)

    if (animate)
      scatter!(ct, [xold[1]],[xold[2]], markercolor=:red, legend=false,
               xlims=xrange, ylims=yrange) 
      quiver!(ct, [xold[1]],[xold[2]], quiver=([α*Δ[1]],[α*Δ[2]]), legend=false,
              xlims=xrange, ylims=yrange)
      frame(anim)
    end

    xchange = norm(x-xold)
    fchange = abs(fnew-fold)
    μiter += 1

    # update μ (the details here could also be improved)    
    foc = ∇ₓL(x,λ)
    if (μiter>10 || (norm(foc)< μ && λ'*c(x)<10*μ)) 
      μ *=  μfactor
      μiter = 0
    end
    
    xold = x
    fold = fnew
    if verbosity>0
      print("Iter $iter: f=$fnew, λ=$λ, c(x)=$(c(x)), μ=$μ, norm(foc)=$(norm(foc))\n")
    end
    iter += 1    
  end
  if (iter >= maxiter)
    info = "Maximum iterations reached"
  elseif (stuck>0)
    info = "Failed to find feasible step for " * string(stuck) * " iterations."
  else
    info = "Convergence."
  end
  return(fold, xold, iter, info, anim) 
end



"""
      sequentialquadratic(f, x0, c; 
                          ∇f = x->ForwardDiff.gradient(f,x),
                          ∇c = x->ForwardDiff.jacobian(c,x),
                          L   = (x,λ)->(f(x) - dot(λ,c(x))),
                          ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
                          ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
                          tol=1e-4, maxiter = 1000,
                          trustradius=1.0, xrange=[-2., 3.],
                          yrange=[-2.,6.], animate=true, verbosity=1)

  
Find the minimum of function `f` subject to `c(x) ≥ 0` by sequential quadratic programming.
  
# Arguments
  
- `f` function to minimizie
- `x0` starting value. Must have c(x0) > 0
- `c` constraint function. Must return an array.
- `∇f = x->ForwardDiff.gradient(f,x)`
- `∇c = x->ForwardDiff.jacobian(c,x)` Jacobian of constraints
- `L   = (x,λ)->(f(x) - dot(λ,c(x)))` Lagrangian
- `∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x)` Derivative of Lagrangian wrt `x`
- `∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x)` Hessian of Lagrangian wrt `x`
- `tol` convergence tolerance
- `maxiter`
- `trustradius` initial trust region radius
- `xrange` range of x-axis for animation
- `yrange` range of y-axis for animation
- `animate` whether to create an animation (if true requires length(x)==2)
- `verbosity` higher values result in more printed output during search. 0 for no output, any number > 0 for some.  
  
# Returns

- `(fmin, xmin, iter, info, animate)` tuple consisting of minimal function
  value, minimizer, number of iterations, and convergence info

"""
function sequentialquadratic(f, x0, c;
                             ∇f = x->ForwardDiff.gradient(f,x),
                             ∇c = x->ForwardDiff.jacobian(c,x),
                             L   = (x,λ)->(f(x) - dot(λ,c(x))),
                             ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
                             ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
                             tol=1e-4, maxiter = 1000,
                             trustradius=1.0,
                             xrange=[-2., 3.],
                             yrange=[-2.,6.], animate=true, verbosity=1)
  fold = f(x0)
  xold = x0
  xchange=Inf
  fchange=Inf
  iter = 0
  μiter = 0
  stuck=0

  animate = animate && length(x0)==2
  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    ct = contour(range(xrange[1],xrange[2], length=100), 
                range(yrange[1],yrange[2], length=100),
                 (x,y) -> log(f([x,y])))
    plot!(ct, xrange, 2.5 .- xrange) # add constraint 
    anim = Animation()
  end
  Dc = ∇c(xold)
  Df = ∇f(xold)
  λ = (Dc*Dc') \ Dc*Df
  foc = ∇ₓL(xold,λ)
  fold  = f(xold)
  negsquared(x) = x < 0 ? x^2 : zero(x)
  merit(x) = f(x) + 1000.0*sum(negsquared.(c(x)))
  while(iter < maxiter && ((xchange>tol) || (fchange>tol) || (stuck>0)
                           || norm(foc)>tol) )
    Df = ∇f(xold)
    Dc = ∇c(xold)
    cx = c(xold)
    H = ∇²ₓL(xold,λ)

    # QP will fail with if H not positive definite
    if !(isposdef(H))
      v,Q = eigen(H)
      Hp = Symmetric(Q*Diagonal(max.(v,one(eltype(v))*1e-8))*Q')
    else
      Hp = H
    end
    # set up and solve our QP
    Δ = Variable(length(xold))
    problem = minimize(Df'*Δ + 0.5*quadform(Δ,Hp), [cx + Dc*Δ >= 0; norm(Δ)<=trustradius])
    solve!(problem, ECOS.Optimizer(verbose=verbosity))
    λ .= problem.constraints[1].dual
    xnew = xold .+ Δ.value

    if (animate)
      scatter!(ct, [xold[1]],[xold[2]], markercolor=:red, legend=false,
               xlims=xrange, ylims=yrange) 
      quiver!(ct, [xold[1]],[xold[2]], quiver=([Δ.value[1]],[Δ.value[2]]), legend=false,
              xlims=xrange, ylims=yrange)
      frame(anim)
    end


    # decide whether to accept new point and whether to adjust trust region
    if (merit(xnew) < merit(xold))
      stuck = 0
      foc = [∇ₓL(xold,λ); λ.*c(xold)]
      if (problem.constraints[2].dual>1e-4) # trust region binding
        trustradius *= 3/2
      end
    else
      stuck += 1
      trustradius *= 2/3
      if (stuck>=20)
        break
      end
    end
    

    xchange = norm(xnew-xold)
    fchange = abs(f(xnew)-f(xold))
    xold = xnew
    
    if verbosity>0
      print("Iter $iter: f=$(f(xold)), λ=$λ, c(x)=$(c(xold)), TR=$trustradius, norm(foc)=$(norm(foc))\n")
    end
    iter += 1    
  end
  if (iter >= maxiter)
    info = "Maximum iterations reached"
  elseif (stuck>0)
    info = "Failed to find feasible step for " * string(stuck) * " iterations."
  else
    info = "Convergence."
  end
  return(f(xold), xold, iter, info, anim) 
end


"""
      slqp(f, x0, c; 
           ∇f = x->ForwardDiff.gradient(f,x),
           ∇c = x->ForwardDiff.jacobian(c,x),
           L   = (x,λ)->(f(x) - dot(λ,c(x))),
           ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
           ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
           tol=1e-4, maxiter = 1000,
           trustradius=1.0, xrange=[-2., 3.],
           yrange=[-2.,6.], animate=true, verbosity=1)

  
Find the minimum of function `f` subject to `c(x) ≥ 0` by sequential
linear quadratic programming. 

See
https://en.wikipedia.org/wiki/Sequential_linear-quadratic_programming
for algorithm information.
  
# Arguments
  
- `f` function to minimizie
- `x0` starting value. Must have c(x0) > 0
- `c` constraint function. Must return an array.
- `∇f = x->ForwardDiff.gradient(f,x)`
- `∇c = x->ForwardDiff.jacobian(c,x)` Jacobian of constraints
- `L   = (x,λ)->(f(x) - dot(λ,c(x)))` Lagrangian
- `∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x)` Derivative of Lagrangian wrt `x`
- `∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x)` Hessian of Lagrangian wrt `x`
- `tol` convergence tolerance
- `maxiter`
- `trustradius` initial trust region radius
- `xrange` range of x-axis for animation
- `yrange` range of y-axis for animation
- `animate` whether to create an animation (if true requires length(x)==2)
- `verbosity` higher values result in more printed output during search. 0 for no output, any number > 0 for some.  
  
# Returns

- `(fmin, xmin, iter, info, animate)` tuple consisting of minimal function
  value, minimizer, number of iterations, and convergence info

"""
function slqp(f, x0, c;
              ∇f = x->ForwardDiff.gradient(f,x),
              ∇c = x->ForwardDiff.jacobian(c,x),
              L   = (x,λ)->(f(x) - dot(λ,c(x))),
              ∇ₓL = (x,λ)->ForwardDiff.gradient(z->L(z,λ), x),
              ∇²ₓL= (x,λ)->ForwardDiff.hessian(z->L(z,λ), x),
              tol=1e-4, maxiter = 1000,
              trustradius=1.0,
              xrange=[-2., 3.],
              yrange=[-2.,6.], animate=true, verbosity=1)
  fold = f(x0)
  xold = x0
  xchange=Inf
  fchange=Inf
  iter = 0
  μiter = 0
  stuck=0
  lptrustradius = trustradius
  
  animate = animate && length(x0)==2
  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    ct = contour(range(xrange[1],xrange[2], length=100), 
                range(yrange[1],yrange[2], length=100),
                 (x,y) -> log(f([x,y])))
    plot!(ct, xrange, 2.5 .- xrange) # add constraint 
    anim = Animation()
  end
  Dc = ∇c(xold)
  Df = ∇f(xold)
  cx = c(xold)
  λ = (Dc*Dc') \ Dc*Df
  foc = [∇ₓL(xold,λ); λ.*cx]
  fold  = f(xold)
  negsquared(x) = x < 0 ? x^2 : zero(x)
  merit(x) = f(x) + sum(negsquared.(c(x)))
  while(iter < maxiter && ((xchange>tol) || (fchange>tol) || (stuck>0)
                           || (norm(foc)>tol)) )
    Df = ∇f(xold)
    Dc = ∇c(xold)
    cx = c(xold)
    H = ∇²ₓL(xold,λ)

    # set up and solve linear program for finding active set (binding
    # constraints)
    Δ = Variable(length(xold))
    problem = minimize(Df'*Δ, [cx + Dc*Δ >= 0,
                               Δ <= lptrustradius,
                               -lptrustradius <= Δ])
    solve!(problem, ECOS.Optimizer(verbose=verbosity))
    while (Int(problem.status)!=1 && lptrustradius <=1e4)
      lptrustradius *= 2
      problem = minimize(Df'*Δ, [cx + Dc*Δ >= 0,
                               Δ <= lptrustradius,
                               -lptrustradius <= Δ])
      solve!(problem, ECOS.Optimizer(verbose=verbosity))
    end

    if isa(problem.constraints[1].dual, AbstractArray)
      active = vec(problem.constraints[1].dual .> tol)
    else
      active = [problem.constraints[1].dual .> tol]
    end
    
    # set up and solve our QP
    Δ = Variable(length(xold))
    if !(isposdef(H))
      Hp = zero(H)
    else
      Hp = H
    end
    if any(active) 
      problem = minimize(Df'*Δ + 0.5*quadform(Δ,Hp), [cx[active] + Dc[active,:]*Δ >= 0, norm(Δ)<=trustradius])
      solve!(problem, ECOS.Optimizer(verbose=verbosity))
      λ[active] .= problem.constraints[1].dual
      λ[.!active] .= zero(eltype(λ))
    else
      problem = minimize(Df'*Δ + 0.5*quadform(Δ,Hp), [norm(Δ)<=trustradius])
      solve!(problem, ECOS.Optimizer(verbose=verbosity))
      λ = zero(cx)
    end
    xnew = xold .+ Δ.value

    if (animate)
      scatter!(ct, [xold[1]],[xold[2]], markercolor=:red, legend=false,
               xlims=xrange, ylims=yrange) 
      quiver!(ct, [xold[1]],[xold[2]], quiver=([Δ.value[1]],[Δ.value[2]]), legend=false,
              xlims=xrange, ylims=yrange)
      frame(anim)
    end


    # decide whether to accept new point and whether to adjust trust region
    if (merit(xnew)  < merit(xold))
      stuck = 0
      foc = [∇ₓL(xold,λ); λ.*c(xold)]
      if (problem.constraints[end].dual>1e-4) # trust region binding
        trustradius *= 3/2
      end
    else
      stuck += 1
      trustradius *= 2/3
      if (stuck>=20)
        break
      end
    end

    xchange = norm(xnew-xold)
    fchange = abs(f(xnew)-f(xold))

    xold = xnew
    
    if verbosity>0
      print("Iter $iter: f=$(f(xold)), λ=$λ, c(x)=$(c(xold)),TR=$trustradius, norm(foc)=$(norm(foc)), $xchange , $fchange, $stuck\n")
    end
    iter += 1    
  end
  if (iter >= maxiter)
    info = "Maximum iterations reached"
  elseif (stuck>0)
    info = "Failed to find feasible step for " * string(stuck) * " iterations."
  else
    info = "Convergence."
  end
  return(f(xold), xold, iter, info, anim) 
end
