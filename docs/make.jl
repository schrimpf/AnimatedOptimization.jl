using Documenter, AnimatedOptimization, DocumenterMarkdown
runweave=false
if runweave
  using Weave

  cd("jmd")
  try
    weave("optimization.jmd",out_path="../build", cache=:user, cache_path="../weavecache",
          doctype="github", mod=Main,
          pandoc_options=["--toc","--toc-depth=2","--filter=pandoc-citeproc"])
    #notebook("optimization.jmd", nbconvert_options="--allow-errors")
  finally 
    cd("..")
  end
  #run(`pandoc jmd/optimization.md -o build/optimization.md -f markdown -t gfm`)
  #cp() # copy gifs to build
  if (isfile("src/temp.md"))
    rm("src/temp.md")
  end
end

makedocs(
  modules=[AnimatedOptimization],
  format=Markdown(),
  clean=false,
  pages=[
    "Home" => "index.md"
  ],
  repo="https://github.com/schrimpf/AnimatedOptimization/blob/{commit}{path}#L{line}",
  sitename="AnimatedOptimization.jl",
  authors="Paul Schrimpf <paul.schrimpf@gmail.com>"
)

run(`mkdocs build`)

deploy=false
if deploy || "deploy" in ARGS

  include("faketravis.jl")

  # deploydocs(deps = nothing, make = nothing,
  #            repo = "github.com/schrimpf/AnimatedOptimization",
  #            target = "build",
  #            branch = "gh-pages",
  #            devbranch = "master"
  #            )

  deploydocs(
    repo = "github.com/schrimpf/AnimatedOptimization",
    deps   = Deps.pip("mkdocs", "pygments", "python-markdown-math", "bootswatch"),
    make   = () -> run(`mkdocs build`),
    target = "site",
    branch = "gh-pages",
    devbranch = "master"
  )

end



## To generate ssh key for Travis
#using DocumenterTools, AnimatedOptimization
#DocumenterTools.genkeys(AnimatedOptimization)
#see https://juliadocs.github.io/Documenter.jl/stable/man/hosting/
