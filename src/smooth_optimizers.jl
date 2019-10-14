"""
     graddescent(f, x0; grad=x->Forwardiff.gradient(f,x),
                 γ0=1.0, ftol = 1e-6,
                 xtol = 1e-4, gtol=1-6, maxiter = 1000, 
                 xrange=[-2., 3.],
                 yrange=[-2.,6.], animate=true)

Find the minimum of function `f` by gradient descent
  
# Arguments
  
- `f` function to minimize
- `x0` starting value
- `grad` function that computes gradient of `f`
- `γ0` initial step size multiplier
- `ftol` tolerance for function value
- `xtol` tolerance for x
- `gtol` tolerance for gradient. Convergence requires meeting all three tolerances.
- `maxiter` maximum iterations
- `xrange` x-axis range for animation
- `yrange` y-axis range for animation
- `animate` whether to create animation

# Returns

- `(fmin, xmin, iter, info, anim)` tuple consisting of minimal function
  value, minimizer, number of iterations, convergence info, and animations
"""
function graddescent(f, x0; grad=x->ForwardDiff.gradient(f,x),
                     γ0=1.0, ftol = 1e-6,
                     xtol = 1e-4, gtol=1-6, maxiter = 1000, 
                     xrange=[-2., 3.],
                     yrange=[-2.,6.], animate=true)
  fold = f(x0)
  xold = x0
  xchange=Inf
  fchange=Inf
  γ = γ0
  iter = 0
  stuck=0
  improve = 0 # we increase γ if 5 steps in a row improve f(x)

  animate = animate && (length(x0)==2)

  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    c = contour(range(xrange[1],xrange[2], length=100),
                range(yrange[1],yrange[2], length=100),
                (x,y) -> log(f([x,y])))
    anim = Animation()
  else
    anim = nothing
  end
  g = grad(xold)
  
  while(iter < maxiter && ((xchange>xtol) || (fchange>ftol) || (stuck>0)
                           || norm(g)>gtol) )
    g = grad(xold)
    x = xold - γ*g
    fnew = f(x)

    if animate
      scatter!(c, [xold[1]],[xold[2]], markercolor=:red, legend=false,
               xlims=xrange, ylims=yrange) 
      quiver!(c, [xold[1]],[xold[2]], quiver=([-γ*g[1]],[-γ*g[2]]), legend=false,
              xlims=xrange, ylims=yrange)
      frame(anim)
    end
    
    if (fnew>=fold)
      γ*=0.5
      improve = 0     
      stuck += 1
      if (stuck>=10)
        break
      end
    else
      stuck = 0
      improve += 1
      if (improve>5)
        γ *= 2.0
        improve=0
      end
      xold = x
      fold = fnew
    end
    xchange = norm(x-xold)
    fchange = abs(fnew-fold)
    iter += 1
  end
  if (iter >= maxiter)
    info = "Maximum iterations reached"
  elseif (stuck>0)
    info = "Failed to improve for " * string(stuck) * " iterations."
  else
    info = "Convergence."
  end
  return(fold, xold, iter, info, anim) 
end

"""
    newton(f, x0; 
           grad=x->ForwardDiff.gradient(f,x),
           hess=x->ForwardDiff.hessian(f,x),
           ftol = 1e-6,
           xtol = 1e-4, gtol=1-6, maxiter = 1000, 
           xrange=[-2., 3.],
           yrange=[-2.,6.], animate=true)

Find the minimum of function `f` by Newton's method.
  
# Arguments
  
- `f` function to minimizie
- `x0` starting value
- `grad` function that returns gradient of `f`
- `hess` function that returns hessian of `f`
- `ftol` tolerance for function value
- `xtol` tolerance for x
- `gtol` tolerance for gradient. Convergence requires meeting all three tolerances.
- `maxiter` maximum iterations
- `xrange` x-axis range for animation
- `yrange` y-axis range for animation
- `animate` whether to create animation

# Returns

- `(fmin, xmin, iter, info, anim)` tuple consisting of minimal function
  value, minimizer, number of iterations, convergence info, and animation

"""
function newton(f, x0;
                grad=x->ForwardDiff.graident(f,x),
                hess=x->ForwardDiff.hession(f,x),
                ftol = 1e-6,
                xtol = 1e-4, gtol=1-6, maxiter = 1000, 
                xrange=[-2., 3.],
                yrange=[-2.,6.], animate=true)
  fold = f(x0)
  xold = x0
  xchange=Inf
  fchange=Inf
  iter = 0
  stuck=0

  animate=animate && length(x0)==2

  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    c = contour(range(xrange[1],xrange[2], length=100),
                range(yrange[1],yrange[2], length=100),
                (x,y) -> log(f([x,y])))
    anim = Animation()
  end
  
  g = grad(xold)
  
  while(iter < maxiter && ((xchange>xtol) || (fchange>ftol) || (stuck>0)
                           || norm(g)>gtol) )
    g = grad(xold)
    H = hess(xold)
    Δx = - inv(H)*g
    x = xold + Δx
    fnew = f(x)

    if animate
      scatter!(c, [xold[1]],[xold[2]], markercolor=:red, legend=false,
               xlims=xrange, ylims=yrange) 
      quiver!(c, [xold[1]],[xold[2]], quiver=([Δx[1]],[Δx[2]]), legend=false,
              xlims=xrange, ylims=yrange)
      frame(anim)
    end
    
    if (fnew>=fold)
      stuck += 1
      if (stuck>=10)
        break
      end
    else
      stuck = 0
      xold = x
      fold = fnew
    end
    xchange = norm(x-xold)
    fchange = abs(fnew-fold)
    iter += 1
  end
  if (iter >= maxiter)
    info = "Maximum iterations reached"
  elseif (stuck>0)
    info = "Failed to improve for " * string(stuck) * " iterations."
  else
    info = "Convergence."
  end
  return(fold, xold, iter, info, anim) 
end
