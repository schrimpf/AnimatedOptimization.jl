using Documenter, AnimatedOptimization, DocumenterMarkdown
runweave=true
if runweave
  using Weave

  cd("build")
  try
    weave("../jmd/optimization.jmd",out_path="../build", cache=:user,
          cache_path="../weavecache", doctype="github", mod=Main,
          args=Dict("md" => "true"))
    #notebook("../jmd/optimization.jmd", out_path="../build",
    #         nbconvert_options="--allow-errors")
    #convert_doc("../jmd/optimization.jmd", "../build/optimization.ipynb")
    #weave("../jmd/optimization.jmd",out_path="../build", cache=:user,
    #      cache_path="../weavecache", doctype="github", mod=Main,
    #      args=Dict("notebook" => "false"))    
  finally 
    cd("..")
  end
  if (isfile("build/temp.md"))
    rm("build/temp.md")
  end
end

makedocs(
  modules=[AnimatedOptimization],
  format=Markdown(),
  clean=false,
  pages=[
    "Home" => "index.md"
  ],
  repo="https://github.com/schrimpf/AnimatedOptimization.jl/blob/{commit}{path}#L{line}",
  sitename="AnimatedOptimization.jl",
  authors="Paul Schrimpf <paul.schrimpf@gmail.com>"
)

run(`mkdocs build`)

deploy=true
if deploy || "deploy" in ARGS
  
  #include("faketravis.jl")

  # deploydocs(deps = nothing, make = nothing,
  #            repo = "github.com/schrimpf/AnimatedOptimization",
  #            target = "build",
  #            branch = "gh-pages",
  #            devbranch = "master"
  #            )

  #deploydocs(
  #  repo = "github.com/schrimpf/AnimatedOptimization",
  #  deps   = Deps.pip("mkdocs", "pygments", "python-markdown-math", "bootswatch"),
  # make   = () -> run(`mkdocs build`),
  # target = "site",
  #  branch = "gh-pages",
  #  devbranch = "master"
  #)
  run(`mkdocs gh-deploy`)
  
end



## To generate ssh key for Travis
#using DocumenterTools, AnimatedOptimization
#DocumenterTools.genkeys(AnimatedOptimization)
#see https://juliadocs.github.io/Documenter.jl/stable/man/hosting/
