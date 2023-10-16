using Random, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 200  #diametro del circulo
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

function barrera(v::Int64, n_pasos::Int64, n_particulas::Int64, radio, Ω=0)
        #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
            x   = fill(NaN,n_pasos,n_particulas)
            y   = similar(x)
            φ   = similar(x)
            forces1 = zeros(n_pasos,n_particulas)
            forces2 = zeros(n_pasos,n_particulas)
            φ[1,:] = 2pi * randn(n_particulas)
            random = sqrt(rand())
            rand_ang = randn(n_particulas)*2pi
            x[1,:], y[1,:] = condicion_inicial(n_particulas,radio,radio_c) 
            @showprogress "Calculando trayectorias " for i in 2:n_pasos
                


                f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], radio)
                forces1[i-1,:] = f_x
                forces2[i-1,:] = f_y 
                x[i-1,:] += f_x*dt
                y[i-1,:] += f_y*dt

                ruidoDtx = sqrtD * randn(n_particulas)
                
                ruidoDty = sqrtD * randn(n_particulas)
                
                ruidoDr  = sqrtT * randn(n_particulas)
                
                φ[i,:] = φ[i-1,:] .+ Ω*dt +  ruidoDr
                
                x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
                
                y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty

                dx = x[i,:]  .- centro_x
                dy = y[i,:]  .- centro_y
                

                #Verificar la posicion de la particula con respecto al centro del circulo
                distancia_c = sqrt.((dx).^2 + (dy).^2)

                # Se define un array con el indice de todas las particulas que están fuera de la barrera en el i-esimo paso temporal
                fuera = distancia_c .>= radio_c

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
                        magnitud_fuerza = lj_fuerza(r,10,2*radio)
                    #Calculamos la componente x e y de la fuerza    
                        f_x = magnitud_fuerza * dx 
                        f_y = magnitud_fuerza * dy
                    # Updateamos el array x e y de las fuerzas
                        fuerza_x[i] += f_x
                        fuerza_y[i] += f_y
                end
            end
        end
    end
    
    return fuerza_x, fuerza_y
end

function lj_fuerza(distancia, epsilon, sigma)
    # Calculate the Lennard-Jones potential energy between two particles
    r_sigma = sigma / distancia

    # Calculate the potential energy
    force = 24 * epsilon * (  ( (sigma^6)/(distancia.^8) ) -2*( (sigma^12)/(distancia.^14) )  )

    return force
end


function condicion_inicial(n_particulas, radio_particula, radio_circulo)
    x_ini = Float64[]  # Array to store x-coordinates
    y_ini = Float64[]  # Array to store y-coordinates

    for _ in 1:n_particulas
        while true
            # Generar condiciones iniciales en coordenadas polares
            angulo = rand() * 2 * π
            r = rand() * (radio_circulo - radio_particula) + radio_particula

            # Convertir a coordenas cartesianas
            x = r * cos(angulo)
            y = r * sin(angulo)

            # Checkear si la posicion de la i-esima particula no se solapa con otra
            overlap = false
            for i in 1:length(x_ini)
                if sqrt((x - x_ini[i])^2 + (y - y_ini[i])^2) < 2 * radio_particula
                    overlapping = true
                    break
                end
            end
            #Si no hay overlap, entonces guarda la posicion en x e y
            if !overlap
                push!(x_ini, x)
                push!(y_ini, y)
                break
            end
        end
    end

    return x_ini, y_ini
end