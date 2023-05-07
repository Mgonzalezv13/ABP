using Random, Serialization, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

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

#Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula ademas del msd y el paso temporal
sqrtD = sqrt(2*Dt*Dr*dt) #esto corresponde a √(2Dt*Dr*dt)

for j in 1:n_particulas
    # Inicializar las variables para cada particula
    x   = zeros(n_trayectorias)
    y   = similar(x)
    φ   = similar(x)
    Δt  = zeros(n_trayectorias)
    msd = similar(Δt)
    x[1]    = 0
    y[1]    = 0
    φ[1]    = 0
    Δt[1]   = 0

    for i in 1:n_trayectorias-1
        ruido  = sqrtD*randn() 
        x[i+1] = x[i] .+ v*cos(φ[i])*dt .+ ruido
        y[i+1] = y[i] .+ v*sin(φ[i])*dt .+  ruido
        φ[i+1] = φ[i] .+ Ω*dt .+  ruido 
        Δt[i+1]  =+ Δt[i] .+ dt
        msd[i+1] =(x[i+1] .- x[1]).^2 + (y[i+1] .- y[1]).^2
    end

    # Save the trajectory to file
    file_path = joinpath(folder_path, "traj_$j.dat")
    writedlm(file_path, [x y])
    file_path = joinpath(folder_path, "msd_$j.dat")
    writedlm(file_path, msd)

    # Plot the MSD and save the figure
    plot(Δt, msd, yaxis=log)
    savefig(joinpath(folder_path, "msd$j.png"))
end









