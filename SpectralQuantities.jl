using Plots,Interpolations,QuadGK,LaTeXStrings,Peaks,ProgressMeter
include("Lx.jl")
##
σ=0.8;
D=σ^2/2
s=1.5
f(x)=s*tanh(x) -x;
L=5;
dx=0.01;
xs=range(-L,stop=L,step=dx);
Us=-cumsum(f.(xs)).*dx # Compute potential function
U=linear_interpolation(xs,Us,extrapolation_bc = 0.0);
ϕ0(x)=exp(-U.(x)./D)
z0=quadgk(ϕ0,-L,L)[1] # Normalization constant
plot(xs,ϕ0.(xs)./z0,xlabel=L"x",ylabel=L"ϕ_0(x)",label=false)
# Testing on ϕ₀:
Nx=1000;
ϕ_q,ψ_δψ=Lx.defineSpectralProblem(f,σ,L,-L,Nx);
ϕt0,qt0,dX,xt=ϕ_q(1,0);
zt0=sum(ϕt0)*dX
p0=plot(xs,ϕ0.(xs)./z0,xlabel=L"x",ylabel=L"ϕ_0(x)",label="Exact")
plot!(xt,real.(ϕt0./zt0),xlabel=L"x",ylabel=L"ϕ_0(x)",label="Numerical",style=:dash)
##
# Find the λs imposing ∂ₓψₙ(±L)=0
function CE(λ)
    ψ,δψ=ψ_δψ(1,λ)
    return real(δψ[end])
end

dλ=0.001
λt=range(-dλ,-6,step=-dλ)
cet=zeros(length(λt));
Threads.@threads for n=1:length(cet)
    cet[n]=abs(CE(λt[n]))
end

λs=λt[argminima(cet)]
p1=scatter(λs,xlabel=L"$n$",ylabel=L"$\lambda_n$",color=:black,label=false)

## Plot some eignfunctions
cs=palette([:red,:blue],4)
p2=plot()

for n=1:4
    ϕt,qt,dX,xt=ϕ_q(1,λs[n]);
    Max=maximum(real.(ϕt))
    plot!(xt,real.(ϕt./Max),xlabel=L"$x$",ylabel=L"$\phi_n(x)$",label="n="*string(n),color=cs[n])
end
plot!()

plot(p0,p1,p2,layout=(1,3),size=(800,300))
savefig("Res.pdf")

