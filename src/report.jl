"Render result tables as a PDF using an external Typst executable."
function write_report(r::AnalysisResult,path;title="BioExplorer community analysis",overwrite=false,compiler="typst")
 endswith(lowercase(path),".pdf")||throw(ArgumentError("Destination must be a PDF"))
 ispath(path)&&!overwrite&&throw(ArgumentError("Refusing to overwrite destination"))
 executable=Sys.which(compiler)
 executable===nothing&&throw(ArgumentError("Install Typst and add its executable to PATH"))
 # User-controlled values are encoded as Typst strings, never as markup.
 literal(x)="\""*replace(string(x),"\\"=>"\\\\","\""=>"\\\"","\n"=>"\\n","\r"=>"\\r","\t"=>"\\t")*"\""
 function section(io,name,table)
  println(io,"== ",name)
  if ncol(table)==0 || nrow(table)==0
   println(io,"No results available.\n");return
  end
  println(io,"#table(columns: ",ncol(table),", inset: 4pt,")
  for name in names(table);println(io,"text(",literal(name),"),");end
  for row in eachrow(table),value in row
   display=ismissing(value) ? "NA" : value isa AbstractFloat ? string(round(value;digits=4)) : string(value)
   println(io,"text(",literal(display),"),")
  end
  println(io,")\n")
 end
 src=sprint() do io
  println(io,"#set page(paper: \"a4\", flipped: true, margin: 1.5cm)")
  println(io,"#set text(size: 9pt)")
  println(io,"= #text(",literal(title),")")
  println(io,"IN DEVELOPMENT: independently validate methods and input sampling design before interpretation.")
  println(io,"Zeros represent unrecorded taxa, not verified absence. Accumulation SD is not a confidence interval.")
  section(io,"Alpha diversity",r.alpha)
  section(io,"Beta dissimilarity",r.beta)
  section(io,"Taxon composition",r.composition)
  section(io,"Site accumulation",r.accumulation)
  r.estimators!==nothing&&section(io,"Abundance estimators",r.estimators)
  section(io,"Input diagnostics",r.diagnostics)
 end
 destination=abspath(path);mkpath(dirname(destination))
 mktempdir() do dir
  source=joinpath(dir,"report.typ");pdf=joinpath(dir,"report.pdf")
  write(source,src);run(`$executable compile $source $pdf`)
  cp(pdf,destination;force=overwrite)
 end
 path
end
