using Random, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 10000 #se bajo el numero de pasos temporales, debido a que con muchos colapsaba el grafico
n_particulas = 100


v  = 0    #velocidad activa de la particula
Dt = 0.22   #Difusion Traslacional
Dr = 0.16   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10e-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)


 for j in 1:n_particulas
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x   = zeros(n_trayectorias)
        y   = similar(x)
        φ   = similar(x)
        x[1]    = 0
        y[1]    = 0
        φ[1]    = 0
        for i in 1:n_trayectorias-1
            
            ruidoDtx  = sqrtD*randn()
            
            ruidoDty  = sqrtD*randn() #Se definen ruidos diferentes para X e Y
            
            ruidoDr  = sqrtT*randn()
            
            φ[i+1] = φ[i] .+ Ω*dt .+  ruidoDr
            
            x[i+1] = x[i] .+ v*cos(φ[i])*dt .+ ruidoDtx
            
            y[i+1] = y[i] .+ v*sin(φ[i])*dt .+  ruidoDty

        end

    #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
    file_path = joinpath(folder_path, "traj_$j.dat")
    writedlm(file_path, [x y])
end









