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
epsilon = 1.0  # Depth of the potential well
sigma = 1.0*2^(1/6)     # Distance at which potential energy is zero
radio_particulas = 1.0

threshold_distance = (2^(1/6) ) * 2 * radio_particulas

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
            rand_ang = randn(n_particulas)*2pi
            x[1,:] = radio*sqrt.(rand(n_particulas)).*cos.(rand_ang)
            y[1,:] = radio*sqrt.(rand(n_particulas)).*sin.(rand_ang)
            @showprogress "Calculando trayectorias " for i in 1:n_pasos-1
                
                ruidoDtx = sqrtD * randn(n_particulas)
                
                ruidoDty = sqrtD * randn(n_particulas)
                
                ruidoDr  = sqrtT * randn(n_particulas)
                
                φ[i+1,:] = φ[i,:] .+ Ω*dt +  ruidoDr
                
                x[i+1,:] = x[i,:] + v*cos.(φ[i,:])*dt + ruidoDtx
                
                y[i+1,:] = y[i,:] + v*sin.(φ[i,:])*dt +  ruidoDty

                dx = x[i+1,:] .- centro_x
                dy = y[i+1,:] .- centro_y
                
               
            #Verificar la posicion de la particula con respecto al centro del circulo
                distancia = sqrt.((x[i+1,:] .- centro_x).^2 + (y[i+1,:] .- centro_y).^2)

                # Se define un array con el indice de todas las particulas que están fuera de la barrera en el i-esimo paso temporal
                fuera = distancia .>= radio

                # Reflect particles outside the barrier
                for j in findall(fuera)
                    #se calcula el vector normal al desplazamiento en x e y
                    normal_x = dx[j] / distancia[j]
                    normal_y = dy[j] / distancia[j]

                    # Reflejar la posicion de la particula
                    reflejo_x = centro_x + normal_x * radio
                    reflejo_y = centro_y + normal_y * radio
                    x[i + 1, j] = reflejo_x
                    y[i + 1, j] = reflejo_y

                    # Refejar el angulo de la particula
                    angulo = atan(dy[j], dx[j])
                    φ[i + 1, j] = angulo + pi / 4
                end
            end
                

    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",y , ',')
    return x, y
end

function lennard_jones_force(distance, epsilon, sigma)
    # Calculate the Lennard-Jones force between two particles
    r_over_sigma_6 = sigma/distance^6
    r_over_sigma_12 = sigma/distance^12

    # Calculate the force magnitude
    force_magnitude =  epsilon * (r_over_sigma_12 - 2*r_over_sigma_6)

    return force_magnitude
end

x, y = barrera(7,10000,500)