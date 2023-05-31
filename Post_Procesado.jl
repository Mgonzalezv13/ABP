using DelimitedFiles, Plots, ProgressMeter, Statistics

folder_path = "/home/mayron/Datos/datos_1"


dt = 10e-3

function msd_prom(v,n_pasos,n_particulas)
    
    msd_sum = zeros(n_pasos,n_particulas)
    Δt = range(dt, length=n_pasos, step=dt)
            

    @showprogress "Calculando..." for i in 1:n_particulas
        msd = readdlm("/home/mayron/Datos/datos_$v/msd_$i.dat")
        msd_sum[:,i] = msd
    end

    msd_promedio = mean(msd_sum, dims=2)
    writedlm("/home/mayron/Datos/datos_$v/msd_promedio$v.csv", msd_promedio, ',')

    plot(Δt, msd_promedio, xaxis=:log, yaxis=:log, xlim=(0.01,1000), ylim=(0.01,10000), size=(2000,8000))
        
end