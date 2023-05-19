using Random, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 10000 #se bajo el numero de pasos temporales, debido a que con muchos colapsaba el grafico
n_particulas = 100

L = 200  #diametro del circulo
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
radio = L/2  # radio del circulo
inicio_gap = 0  # inicio del agujero en la  barrera circular (en radianes)
fin_gap = 0 # fin del agujero en la barrera circular (en radianes)

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

            # Angulo de la particula con respecto al centro del circulo
            dx = x[i+1] - centro_x
            dy = y[i+1] - centro_y
            angulo = atan(dy, dx)

             Se ajusta el angulo entre -pi y pi
            angulo = angulo % (2*pi)
            if angulo > pi
                angulo -= 2*pi
            end

             ver si la particula esta dentro de la zona del gap
            if angulo >= inicio_gap && angulo <= fin_gap
                continue  # Skip reflection if the particle is within the gap region
            end
                # Mientras la particula este dentro del circulo y fuera del angulo del gap
            distancia = sqrt(dx^2 + dy^2)
                
            if distancia >= radio
                    # Reflexion de la particula en la barrera
                    normal_x = dx / distancia  # componente x del vector normal
                    normal_y = dy / distancia  # componente y del vector normal
                    
                    # Actualizacion de la posicion despues de la reflexion
                    reflejo_x = centro_x + normal_x * radio
                    reflejo_y = centro_y + normal_y * radio
                    delta_x = reflejo_x - x[i+1]
                    delta_y = reflejo_y - y[i+1]
                    x[i+1] = reflejo_x + delta_x
                    y[i+1] = reflejo_y + delta_y
                    
                    # Angulo de reflexion
                    angle = atan(delta_y, delta_x)
                    φ[i+1] = angulo + pi + (angulo - φ[i])
            end
        end

    #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
    file_path = joinpath(folder_path, "traj_$j.dat")
    writedlm(file_path, [x y])
end









