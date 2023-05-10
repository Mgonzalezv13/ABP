using Random, Serialization, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 1000000
n_particulas = 500

#Aca se definen las constantes 
L  = 40
v  = 0.0
Dt = 0.22
Dr = 0.16
Ω  = 0.0
dt = 10e-3


sqrtD = sqrt(2*Dt*Dr*dt) #esto corresponde a √(2Dt*Dr*dt)

for j in 1:n_particulas
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
    x   = zeros(n_trayectorias)
    y   = similar(x)
    φ   = similar(x)
    msd = similar(x)
    x[1]    = 0.25*randn()
    y[1]    = 0.3*randn()
    φ[1]    = 0
    for i in 1:n_trayectorias-1
        ruido  = sqrtD*randn() 
        x[i+1] = x[i] .+ v*cos(φ[i])*dt .+ ruido
        y[i+1] = y[i] .+ v*sin(φ[i])*dt .+  ruido
        φ[i+1] = φ[i] .+ Ω*dt .+  ruido 
        msd[i] =(x[i+1] .- x[1]).^2 + (y[i+1] .- y[1]).^2
        msd[1000000] = (x[1000000] .- x[1]).^2 + (y[1000000] .- y[1]).^2
    end

    #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
    #file_path = joinpath(folder_path, "traj_$j.dat")
    #writedlm(file_path, [x y])
    file_path = joinpath(folder_path, "msd_$j.dat")
    writedlm(file_path, msd)
end









