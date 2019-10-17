using Documenter, AnimatedOptimization
runweave=false
if runweave
  using Weave
  cd("src")
  try
    weave("optimization.jmd", cache=:user,
          doctype="pandoc2html", mod=Main,
          pandoc_options=["--toc","--toc-depth=2","--filter=pandoc-citeproc"])
    notebook("optimization.jmd", nbconvert_options="--allow-errors")
  finally 
    cd("..")
  end
  
  if (isfile("src/temp.md"))
    rm("src/temp.md")
  end
end

makedocs(;
    modules=[AnimatedOptimization],
    format=Documenter.HTML(),
    pages=[
      "Home" => "index.md",
      "Animations" => "redirect.md"
    ],
    repo="https://github.com/schrimpf/AnimatedOptimization.jl/blob/{commit}{path}#L{line}",
    sitename="AnimatedOptimization.jl",
    authors="Paul Schrimpf <paul.schrimpf@gmail.com>"
)

deploydocs(;
    repo="github.com/schrimpf/AnimatedOptimization.jl",
)


## To generate ssh key for Travis
#using DocumenterTools, AnimatedOptimization
#Travis.genkeys(AnimatedOptimization)
#see https://juliadocs.github.io/Documenter.jl/stable/man/hosting/
