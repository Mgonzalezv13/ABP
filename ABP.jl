using Distributions
using Plots
import Random

n_trayectorias = 1000000
n_particulas = 10

L  = 40
v  = fill(10.0, n_trayectorias)
Dt = fill(0.08, n_trayectorias)
Dr = fill(0.16, n_trayectorias)
Ω  = zeros(n_trayectorias)
dt = 10e-3

x   = zeros(n_trayectorias, n_particulas)
y   = similar(x)
φ   = similar(x)
Δt  = dt*(1:n_trayectorias)
collect(Δt)
msd = similar(x)
x[1, :] .= 0.0
y[1, :] .= 0.0
φ[1, :] .= 0.0


for j in 1:n_particulas
    for i in 1:n_trayectorias-1
        x[i+1, j] = x[i, j] + v[i] * cos(φ[i, j]) * dt + sqrt.(2 * Dt[i] * dt) .* randn()
        y[i+1, j] = y[i, j] + v[i] * sin(φ[i, j]) * dt + sqrt.(2 * Dt[i] * dt) .* randn()
        φ[i+1, j] = φ[i, j] + Ω[i] * dt + sqrt.(2 * Dr[i] * dt) .* randn()
    end
    msd[:, j] = (x[:, j] .- x[1, j]).^2 + (y[:, j] .- y[1, j]).^2
    plot(Δt, msd[:, j], yaxis = log)
    savefig("msd$j.svg")
end
