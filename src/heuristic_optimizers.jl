"""
    minrandomsearch(f, x0, npoints; var0=1.0, ftol = 1e-6,
                         vtol = 1e-4, maxiter = 1000,
                         vshrink=0.9, xrange=[-2., 3.],
                         yrange=[-2.,6.])

Find the minimum of function `f` by random search.

Creates an animation illustrating search progress.

# Arguments

- `f` function to minimize
- `x0` starting value
- `npoints` number of points in cloud
- `var0` initial variance of points
- `ftol` convergence tolerance for function value. Search terminates if both function change is less than ftol and variance is less than vtol.
- `vtol` convergence tolerance for variance. Search terminates if both function change is less than ftol and variance is less than vtol.
- `maxiter` maximum number of iterations
- `vshrink` after every 100 iterations with no function improvement, the variance is reduced by this factor
- `xrange` range of x-axis in animation
- `yrange` range of y-axis in animation
- `animate` whether to create animation

# Returns

- `(fmin, xmin, iter, info, anim)` tuple consisting of minimal function
  value, minimizer, number of iterations, convergence info, and an animation
"""
function minrandomsearch(f, x0, npoints; var0=1.0, ftol = 1e-6,
                         vtol = 1e-4, maxiter = 1000,
                         vshrink=0.9, xrange=[-2., 3.],
                         yrange=[-2.,6.], animate=true)
  var = var0     # current variance for search
  oldbest = Inf  # smallest function value
  xbest = x0     # x with smallest function vale
  newbest = f(xbest)
  iter = 0       # number of iterations
  noimprove = 0  # number of iterations with no improvement

  animate = (animate && length(x0)==2)

  if animate
    # make a contour plot of the function we're minimizing. This is for
    # illustrating; you wouldn't have this normally
    x = range(xrange[1],xrange[2], length=100)
    y = range(yrange[1],yrange[2], length=100)
    c = contour(x,y,(x,y) -> log(f([x,y])))
    anim = Animation()
  else
    anim = nothing
  end

  while ((oldbest - newbest > ftol || var > vtol) && iter<=maxiter)
    oldbest = newbest
    x = rand(MvNormal(xbest, var),npoints)

    if animate
      # plot the search points
      p = deepcopy(c)
      scatter!(p, x[1,:], x[2,:], markercolor=:black, markeralpha=0.5, legend=false, xlims=xrange, ylims=yrange)
    end

    fval = mapslices(f,x, dims=[1])
    (newbest, i) = findmin(fval)
    if (newbest > oldbest)
      noimprove+=1
      newbest=oldbest
    else
      xbest = x[:,i[2]]
      noimprove = 0
    end

    if animate
      # plot the best point so far
      scatter!(p, [xbest[1]],[xbest[2]], markercolor=:red, legend=false)
    end

    if (noimprove > 10) # shrink var
      var *= vshrink
    end

    if animate
      frame(anim) # add frame to animation
    end

    iter += 1
  end
  if (iter>maxiter)
    info = "Maximum iterations reached"
  else
    info = "Convergence."
  end
  return(minimum = newbest, minimizer = xbest, iters = iter, info = info, anim = anim)
end
