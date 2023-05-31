using Random, DelimitedFiles, Distributions, ProgressMeter




Dt = 0.22   #Difusion Traslacional
Dr = 0.16   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10e-3  #Paso temporal

sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)

 function posicion(v,Ω, n_pasos, n_particulas)
    folder_path = "/home/mayron/Datos/datos_$v"
    mkdir(folder_path)
    @showprogress "Calculando las trayectorias..." for j in 1:n_particulas
            #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
                x   = zeros(n_pasos)
                y   = similar(x)
                φ   = similar(x)
                msd = similar(x)
                x[1]    = randn()
                y[1]    = randn()
                φ[1]    = rand(Uniform(0, 2π))
                        for i in 1:n_pasos-1
                            
                            ruidoDtx  = sqrtD*randn()
                            
                            ruidoDty  = sqrtD*randn() #Se definen ruidos diferentes para X e Y
                            
                            ruidoDr  = sqrtT*randn()
                            
                            φ[i+1] = φ[i] + Ω*dt +  ruidoDr
                            
                            x[i+1] = x[i] + v*cos(φ[i])*dt + ruidoDtx
                            
                            y[i+1] = y[i] + v*sin(φ[i])*dt +  ruidoDty

                            msd[i+1] = ( x[i+1]- x[1] )^2 + ( y[i+1] - y[1] )^2
                        end
                        msd[1]  = ( x[2]- x[1] )^2 + ( y[2] - y[1] )^2

            #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
            file_path = joinpath(folder_path, "traj_$j.dat")
            writedlm(file_path, [x y])
            file_path = joinpath(folder_path, "msd_$j.dat")
            writedlm(file_path, msd)
        end
end