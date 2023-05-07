using DelimitedFiles, Plots

n_trayectorias=1000000
n_particulas=500
dt = 10e-3
msd_sum = zeros(n_trayectorias)
Δt      = similar(msd_sum)
Δt[1]   = dt


for i in 1:n_trayectorias-1
    
    Δt[i+1] = Δt[i] .+ dt
    
end

for i in 1:n_particulas
    
    msd = readdlm("msd_$i.dat")

    for j in 1:n_trayectorias-1
        msd_sum[j] += msd[j]
    end

end
        
    msd_mean = msd_sum ./ n_particulas  # Divide by n_particulas to get the mean MSD
    plot(Δt,msd_mean,xaxis=log, yaxis=log)

