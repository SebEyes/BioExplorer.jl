"Validated site-by-taxon matrix with explicit measurement semantics."
struct CommunityMatrix{T<:Real,M<:AbstractMatrix{T}}
 sites::Vector{String}
 taxa::Vector{String}
 values::M
 semantics::Symbol
 metadata::Dict{String,Any}
 function CommunityMatrix(sites,taxa,values::AbstractMatrix{T};semantics=:abundance,metadata=Dict{String,Any}()) where T<:Real
  semantics in (:incidence,:abundance,:records) || throw(ArgumentError("Invalid semantics"))
  size(values)==(length(sites),length(taxa)) || throw(DimensionMismatch("Labels must match matrix dimensions"))
  ss,tt=String.(sites),String.(taxa)
  length(unique(ss))==length(ss) && length(unique(tt))==length(tt) || throw(ArgumentError("Labels must be unique"))
  all(isfinite,values) && all(>=(0),values) || throw(ArgumentError("Values must be finite and non-negative"))
  semantics==:incidence && !all(x->x==0||x==1,values) && throw(ArgumentError("Incidence values must be binary"))
  new{T,typeof(copy(values))}(ss,tt,copy(values),semantics,Dict{String,Any}(string(k)=>v for (k,v) in metadata))
 end
end
_blank(x)=ismissing(x)||isempty(strip(string(x)))
"Return row-level diagnostics for Darwin Core-style tabular records."
function validate_records(table;event_col=:eventID,taxon_col=:scientificName,status_col=:occurrenceStatus,count_col=nothing)
 d=DataFrame(table); ns=Set(Symbol.(names(d)))
 out=DataFrame(row=Union{Missing,Int}[],column=String[],code=String[],severity=String[],message=String[])
 for c in (event_col,taxon_col)
  c in ns || push!(out,(missing,string(c),"missing_column","error","Required identifier column is absent."))
 end
 if count_col!==nothing && !(count_col in ns)
  push!(out,(missing,string(count_col),"missing_column","error","Requested count column is absent."))
 end
 for i in 1:nrow(d)
  for c in (event_col,taxon_col)
   c in ns && _blank(d[i,c]) && push!(out,(i,string(c),"missing_identifier","error","Missing identifiers are excluded, never grouped."))
  end
  if status_col in ns
   st=d[i,status_col]
   if _blank(st)
    push!(out,(i,string(status_col),"missing_status","warning","Missing status retained; absence is not inferred."))
   elseif lowercase(strip(string(st))) ∉ ("present","absent")
    push!(out,(i,string(status_col),"unknown_status","warning","Unknown status is excluded."))
   end
  end
  if count_col!==nothing && count_col in ns
   x=d[i,count_col]; q=ismissing(x) ? nothing : tryparse(Float64,string(x))
   if q===nothing || !isfinite(q) || q<0
    push!(out,(i,string(count_col),"invalid_count","error","Count must be finite and non-negative; no imputation."))
   end
  end
 end
 out
end
"Build communities and return (community, diagnostics, exclusions). Zeros mean unrecorded, not verified absence."
function from_records(table;event_col=:eventID,taxon_col=:scientificName,status_col=:occurrenceStatus,count_col=nothing,semantics=:incidence)
 semantics in (:incidence,:abundance,:records)||throw(ArgumentError("Invalid semantics"))
 d=DataFrame(table); ns=Set(Symbol.(names(d)))
 diag=validate_records(d;event_col,taxon_col,status_col,count_col)
 all(c->c in ns,(event_col,taxon_col))||throw(ArgumentError("Required identifier columns absent"))
 semantics==:abundance && (count_col===nothing || !(count_col in ns)) && throw(ArgumentError("Abundance requires an explicit count column"))
 sites=String[]; taxa=String[]; vals=Dict{Tuple{String,String},Float64}()
 exc=DataFrame(source_row=Int[],reason=String[])
 for i in 1:nrow(d)
  ev,tx=d[i,event_col],d[i,taxon_col]
  if _blank(ev)||_blank(tx); push!(exc,(i,"missing identifier"));continue;end
  st=status_col in ns ? d[i,status_col] : missing
  if !_blank(st)
   status=lowercase(strip(string(st)))
   if status=="absent";push!(exc,(i,"explicit absence"));continue;end
   if status!="present";push!(exc,(i,"unknown status"));continue;end
  end
  v=1.0
  if semantics==:abundance
   raw=d[i,count_col]; q=ismissing(raw) ? nothing : tryparse(Float64,string(raw))
   if q===nothing||!isfinite(q)||q<0;push!(exc,(i,"invalid or missing count"));continue;end
   v=q
  end
  e,t=strip(string(ev)),strip(string(tx))
  e in sites||push!(sites,e);t in taxa||push!(taxa,t)
  key=(e,t);vals[key]=get(vals,key,0.0)+v
 end
 m=zeros(length(sites),length(taxa));si=Dict(x=>i for(i,x) in enumerate(sites));ti=Dict(x=>i for(i,x) in enumerate(taxa))
 for ((e,t),v) in vals;m[si[e],ti[t]]=semantics==:incidence ? 1.0 : v;end
 meta=Dict("source_format"=>"Darwin Core-style table","zero_interpretation"=>"unrecorded, not verified absence","input_rows"=>nrow(d),"excluded_rows"=>nrow(exc))
 CommunityMatrix(sites,taxa,m;semantics,metadata=meta),diag,exc
end
"Read a CSV or TSV file; full DwC-A archives are not supported."
function from_table(path;delim=nothing,kwargs...)
 isfile(path)||throw(ArgumentError("File not found"))
 sep=delim===nothing ? (endswith(lowercase(path),".tsv") ? '\t' : ',') : delim
 from_records(CSV.read(path,DataFrame;delim=sep);kwargs...)
end
