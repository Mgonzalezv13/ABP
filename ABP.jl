using Random, Plots, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 50  #diametro del circulo
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
radio = L/2  # radio del circulo
inicio_gap = 0
fin_gap =   pi/6
radio1 = L/2 + 1
Dt = 0.0   #Difusion Traslacional
Dr = 0.3   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)
uniform_dist = Uniform(0, 2π)
epsilon = 0.01  # Depth of the potential well
sigma = 0.5     # Distance at which potential energy is zero
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
            forces1 = zeros(n_pasos,n_particulas)
            forces2 = zeros(n_pasos,n_particulas)
            φ[1,:] = 2pi * randn(n_particulas)
            random = sqrt(rand())
            rand_ang = randn(n_particulas)*2pi
            x[1,:] = radio*sqrt.(rand(n_particulas)).*cos.(rand_ang)
            y[1,:] = radio*sqrt.(rand(n_particulas)).*sin.(rand_ang)
            @showprogress "Calculando trayectorias " for i in 2:n_pasos
                


                f_x, f_y = correccion_disco_solido(x[i-1,:], y[i-1,:], 1,100)
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
                distancia = sqrt.((dx).^2 + (dy).^2)

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



function correccion_disco_solido(posicion_x, posicion_y, radio,potencial)
    n_particulas = length(posicion_x)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
    for i in 1:n_particulas
        for j in 1:n_particulas
            if i != j
                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                if r < 2 * radio
                    #Potencial de interaccion de disco solido
                        force_magnitude = -potencial * (2*radio - r)
                    #Calculamos la componente x e y de la fuerza    
                        f_x = force_magnitude * (dx / r)
                        f_y = force_magnitude * (dy / r)
                    # Updateamos el array x e y de las fuerzas
                        fuerza_x[i] += f_x
                        fuerza_y[i] += f_y
                end
            end
        end
    end
    
    return fuerza_x, fuerza_y
end


#x, y = barrera(7,1000,3)