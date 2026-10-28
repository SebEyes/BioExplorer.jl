"Per-site richness, Shannon, Simpson concentration, Hill numbers and Pielou evenness."
function alpha_diversity(c::CommunityMatrix)
 d=DataFrame(site=String[],richness=Int[],shannon=Float64[],simpson=Float64[],hill_q1=Float64[],hill_q2=Float64[],pielou=Union{Missing,Float64}[])
 for (i,s) in enumerate(c.sites)
  x=Float64.(c.values[i,:]);n=sum(x);p=n==0 ? zeros(length(x)) : x./n
  r=count(>(0),x);h=-sum(q>0 ? q*log(q) : 0.0 for q in p);si=sum(p.^2)
  push!(d,(s,r,h,si,r==0 ? 0.0 : exp(h),si==0 ? 0.0 : 1/si,r<=1 ? missing : h/log(r)))
 end
 d
end
"Long pairwise table; zeros denote unrecorded taxa, not confirmed absences."
function beta_diversity(c::CommunityMatrix)
 d=DataFrame(site1=String[],site2=String[],jaccard=Float64[],sorensen=Float64[],bray_curtis=Float64[])
 for i in 1:(length(c.sites)-1),j in (i+1):length(c.sites)
  a=Float64.(c.values[i,:]);b=Float64.(c.values[j,:])
  u=count(k->a[k]>0||b[k]>0,eachindex(a));sh=count(k->a[k]>0&&b[k]>0,eachindex(a));total=sum(a)+sum(b)
  push!(d,(c.sites[i],c.sites[j],u==0 ? 0.0 : 1-sh/u,u+sh==0 ? 0.0 : 1-2sh/(u+sh),total==0 ? 0.0 : sum(abs.(a.-b))/total))
 end
 d
end
"Rank positive taxon values within each site."
function rank_abundance(c::CommunityMatrix)
 d=DataFrame(site=String[],rank=Int[],taxon=String[],value=Float64[],proportion=Float64[])
 for (i,s) in enumerate(c.sites)
  x=c.values[i,:];n=sum(x);n==0&&continue
  for (r,j) in enumerate(sortperm(x;rev=true))
   x[j]>0&&push!(d,(s,r,c.taxa[j],Float64(x[j]),Float64(x[j]/n)))
  end
 end
 d
end
"Randomized site accumulation without input mutation. SD is not a confidence interval."
function accumulation_curve(c::CommunityMatrix;iterations=100,rng=Random.default_rng())
 iterations isa Integer && iterations>0||throw(ArgumentError("iterations must be a positive integer"))
 n=length(c.sites);x=zeros(iterations,n)
 for r in 1:iterations
  p=falses(length(c.taxa))
  for (k,i) in enumerate(randperm(rng,n));p .|= c.values[i,:].>0;x[r,k]=count(p);end
 end
 DataFrame(sites_added=collect(1:n),mean_richness=n==0 ? Float64[] : vec(mean(x,dims=1)),sd_richness=n==0 ? Float64[] : iterations==1 ? zeros(n) : vec(std(x,dims=1)))
end
"Abundance-only Chao1 and Good's singleton sample coverage; requires integer counts."
function richness_estimators(c::CommunityMatrix)
 c.semantics==:abundance||throw(ArgumentError("Chao1 requires abundance semantics"))
 all(isinteger,c.values)||throw(ArgumentError("Chao1 requires integer individual counts"))
 d=DataFrame(site=String[],observed_richness=Int[],singletons=Int[],doubletons=Int[],chao1=Float64[],sample_coverage=Union{Missing,Float64}[])
 for (i,s) in enumerate(c.sites)
  x=Float64.(c.values[i,:]);r=count(>(0),x);f1=count(==(1),x);f2=count(==(2),x)
  push!(d,(s,r,f1,f2,r+(f2==0 ? f1*(f1-1)/2 : f1^2/(2f2)),sum(x)==0 ? missing : 1-f1/sum(x)))
 end
 d
end
sampling_coverage(c::CommunityMatrix)=richness_estimators(c)[!,[:site,:sample_coverage]]
struct AnalysisResult
 alpha::DataFrame
 beta::DataFrame
 composition::DataFrame
 accumulation::DataFrame
 estimators::Union{Nothing,DataFrame}
 diagnostics::DataFrame
end
"Run supported analyses and retain reusable result tables."
analyze(c::CommunityMatrix;iterations=100,rng=Random.default_rng(),diagnostics=DataFrame())=AnalysisResult(alpha_diversity(c),beta_diversity(c),rank_abundance(c),accumulation_curve(c;iterations,rng),c.semantics==:abundance ? richness_estimators(c) : nothing,DataFrame(diagnostics))
"Export result tables, checking all destinations before writing."
function export_tables(r::AnalysisResult,dir;overwrite=false)
 ts=Dict("alpha.csv"=>r.alpha,"beta.csv"=>r.beta,"composition.csv"=>r.composition,"accumulation.csv"=>r.accumulation,"diagnostics.csv"=>r.diagnostics)
 r.estimators!==nothing&&(ts["estimators.csv"]=r.estimators)
 paths=[joinpath(dir,n) for n in keys(ts)]
 !overwrite&&any(ispath,paths)&&throw(ArgumentError("Refusing to overwrite existing tables"))
 mkpath(dir)
 for (n,d) in ts;CSV.write(joinpath(dir,n),d);end
 paths
end
