using DelimitedFiles, Plots

folder_path = "/home/mayron/ABP/"

n_trayectorias = 1000000
n_particulas = 500
dt = 10e-3

msd_sum = zeros(n_trayectorias)
deltat  = dt*(1:n_trayectorias)
Δt      = collect(deltat)

for i in 1:n_particulas
    msd = readdlm("/home/mayron/ABP/msd_$i.dat")

    @. msd_sum += msd  # Use broadcasting for vectorized addition
end
msd_mean = msd_sum / n_particulas
plot(Δt, msd_mean, xaxis=log, yaxis=log)