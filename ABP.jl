using Random, Plots, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 40  #diametro del circulo
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
radio = L/2  # radio del circulo
inicio_gap = 0
fin_gap =   pi/6
radio1 = L/2 + 1
Dt = 0.22   #Difusion Traslacional
Dr = 0.16   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)
uniform_dist = Uniform(0, 2π)

function barrera(v, n_pasos, n_particulas, Ω=0)
    #pos_x = fill(NaN,n_pasos,n_particulas)
    #pos_y = similar(pos_x)
   # @showprogress "Calculando trayectorias " for j in 1:n_particulas
        #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
            x   = fill(NaN,n_pasos,n_particulas)
            y   = similar(x)
            φ   = similar(x)
            φ[1,:] = 2pi * randn(n_particulas)
            random = sqrt(rand())
            rand_ang = randn()*2pi
            x[1,:] = (rand(n_particulas))
            y[1,:] = (rand(n_particulas))
            @showprogress "Calculando trayectorias " for i in 1:n_pasos-1
                
                ruidoDtx = sqrtD * randn(n_particulas)
                
                ruidoDty = sqrtD * randn(n_particulas)
                
                ruidoDr  = sqrtT * randn(n_particulas)
                
                φ[i+1,:] = φ[i,:] .+ Ω*dt +  ruidoDr
                
                x[i+1,:] = x[i,:] + v*cos.(φ[i,:])*dt + ruidoDtx
                
                y[i+1,:] = y[i,:] + v*sin.(φ[i,:])*dt +  ruidoDty

             #   dx = x[i+1] - centro_x
              #  dy = y[i+1] - centro_y
                
               
            #Verificar la posicion de la particula con respecto al centro del circulo
               # distancia = sqrt((x[i+1] - centro_x)^2 + (y[i+1] - centro_y)^2)

            #Verificar si la particula esta fuera de la barrera
                #if distancia >= radio
                    # Reflect the particle off the circular wall
                 #   normal_x = dx / distancia  # x-component of outward normal vector
                  #  normal_y = dy / distancia  # y-component of outward normal vector
                    
                    # Reflect the particle's position
                   # reflejo_x =  normal_x * radio
                   # reflejo_y =  normal_y * radio
                   # delta_x   = reflejo_x - x[i+1]
                   # delta_y   = reflejo_y - y[i+1]
                   # x[i+1]    = reflejo_x + delta_x
                   # y[i+1]    = reflejo_y + delta_y
                    
                    # Reflect the angle
                    #angulo = atan(delta_y, delta_x)
                    #φ[i+1] = angulo + pi/4 
                #end
           # pos_x[:,j] = x
           # pos_y[:,j] = y
            end
                

    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",y , ',')
end