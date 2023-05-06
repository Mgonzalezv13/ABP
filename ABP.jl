using Distributions, LinearAlgebra, GLM, DataFrames, Printf, CSV
using Plots, Plots.PlotMeasures, StatsPlots
import Random


#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 1000000
n_particulas = 10

#Aca se definen las constantes 
L  = 40
v  = 10.0
Dt = 0.08
Dr = 0.16
Ω  = 0.0
dt = 10e-3

#Aca se definen vectores "vacios
x   = zeros(n_trayectorias)
y   = similar(x)
φ   = similar(x)
Δt  = zeros(n_trayectorias)
msd = similar(Δt)
x[1]    = 0
y[1]    = 0
φ[1]    = 0
Δt[1]   = 0

for j in 1:n_particulas
    for i in 1:n_trayectorias-1
            x[i+1] = x[i] .+ v*cos(φ[i])*dt .+ sqrt(2*Dt*Dr*dt)*randn()
            y[i+1] = y[i] .+ v*sin(φ[i])*dt .+ sqrt(2*Dt*Dr*dt)*randn()
            φ[i+1] = φ[i] .+ Ω*dt .+ sqrt(2*Dt*Dr*dt)*randn()
            Δt[i+1]  =+ Δt[i] .+ dt
            msd[i+1] =(x[i] .- x[1]).^2 + (y[i] .- y[1]).^2
    end
    plot(Δt,msd, yaxis=log)
    savefig("msd$j.png")
end














