using DelimitedFiles, Plots


function msd_mean(n_trayectorias,n_particulas)
    
    msd_sum[1] = readdlm("msd_1.dat")

    for i in 1:n_particulas
        msd = readdlm("msd_$i.dat")

        for j in 1:n_trayectorias
            msd_sum[j] += msd[j]
        end

    end
        
    msd_mean = msd_sum ./ n_particulas  # Divide by n_particulas to get the mean MSD
    plot(Δt,msd_mean, yaxis=log)

end