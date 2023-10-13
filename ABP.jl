using Random, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 50  #diametro del circulo
radio_c = L/2
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
inicio_gap = 0
fin_gap =   pi/6
Dt = 0   #Difusion Traslacional
Dr = 0.3   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)
uniform_dist = Uniform(0, 2π)

function barrera(v, n_pasos, n_particulas, radio, Ω=0)
    #pos_x = fill(NaN,n_pasos,n_particulas)
    #pos_y = similar(pos_x)
   # @showprogress "Calculando trayectorias " for j in 1:n_particulas
        #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
            x   = fill(NaN,n_pasos,n_particulas)
            y   = similar(x)
            φ   = similar(x)
            forces1 = zeros(n_pasos,n_particulas)
            forces2 = zeros(n_pasos,n_particulas)
            φ[1,:] = 2pi * randn(n_particulas)
            random = sqrt(rand())
            rand_ang = randn(n_particulas)*2pi
            x[1,:] = radio*sqrt.(rand(n_particulas)).*cos.(rand_ang)
            y[1,:] = radio*sqrt.(rand(n_particulas)).*sin.(rand_ang)
            @showprogress "Calculando trayectorias " for i in 2:n_pasos
                


                f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], radio)
                forces1[i,:] = f_x
                forces2[i,:] = f_y 
                x[i-1,:] += f_x*dt
                y[i-1,:] += f_y*dt

                ruidoDtx = sqrtD * randn(n_particulas)
                
                ruidoDty = sqrtD * randn(n_particulas)
                
                ruidoDr  = sqrtT * randn(n_particulas)
                
                φ[i,:] = φ[i-1,:] .+ Ω*dt +  ruidoDr
                
                x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
                
                y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty

                dx = x[i,:] .- centro_x
                dy = y[i,:] .- centro_y
                

            #Verificar la posicion de la particula con respecto al centro del circulo
                distancia_c = sqrt.((dx).^2 + (dy).^2)

                # Se define un array con el indice de todas las particulas que están fuera de la barrera en el i-esimo paso temporal
                fuera = distancia_c .>= radio

                # Reflect particles outside the barrier
                for j in findall(fuera)
                    #se calcula el vector normal al desplazamiento en x e y
                    normal_x = dx[j] / distancia_c[j]
                    normal_y = dy[j] / distancia_c[j]

                    # Reflejar la posicion de la particula
                    reflejo_x = centro_x + normal_x * radio_c
                    reflejo_y = centro_y + normal_y * radio_c
                    x[i, j] = reflejo_x
                    y[i, j] = reflejo_y

                    # Refejar el angulo de la particula
                    angulo = atan(dy[j], dx[j])
                    φ[i, j] = angulo + pi / 4
                end
            end
                

    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",y , ',')
    writedlm("/home/mayron/Datos/barrera_$v/f_x_v=00$v.csv",forces1 , ',')
    writedlm("/home/mayron/Datos/barrera_$v/f_y_v=00$v.csv",forces2 , ',')
    return x, y
end



function correccion_lj(posicion_x, posicion_y, radio)
    n_particulas = length(posicion_x)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
    for i in 1:n_particulas
        for j in 1:n_particulas
            if i != j
                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                if r <= 2*radio * 2^(1/6)
                    #Potencial de interaccion de disco solido
                        force_magnitude = lj_force(r,100,2*radio)
                    #Calculamos la componente x e y de la fuerza    
                        f_x = force_magnitude * dx 
                        f_y = force_magnitude * dy
                    # Updateamos el array x e y de las fuerzas
                        fuerza_x[i] += f_x
                        fuerza_y[i] += f_y
                end
            end
        end
    end
    
    return fuerza_x, fuerza_y
end

function lj_force(distance, epsilon, sigma)
    # Calculate the Lennard-Jones potential energy between two particles
    r_sigma = sigma / distance

    # Calculate the potential energy
    force = 24 * epsilon * (  ( (sigma^6)/(distance.^8) ) -2*( (sigma^12)/(distance.^14) )  )

    return force
end

function lj_potential(distance, epsilon, sigma)
    # Calculate the Lennard-Jones potential energy between two particles
    r_sigma = sigma / distance

    # Calculate the potential energy
    potential = 4 * epsilon * ((r_sigma.^12) - (r_sigma.^7))

    return potential
end

#x, y = barrera(1,1000,1,1)