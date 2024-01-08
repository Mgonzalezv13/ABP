using Random, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 50  #diametro del circulo
radio_c = L/2
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
#inicio_gap = 0
#fin_gap =   pi/6
Dt = 0   #Difusion Traslacional
Dr = 0.3   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-4  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)
uniform_dist = Uniform(0, 2π)

function barrera(v::Int64, n_pasos::Int64, n_particulas::Int64, radio, Ω=0)
        #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
            Ro = 4*radio
            x   = zeros(n_pasos,n_particulas)
            y   = similar(x)
            φ   = similar(x)
            forces1 = zeros(n_pasos,n_particulas)
            forces2 = zeros(n_pasos,n_particulas)
            φ[1,:] = rand(-pi:pi,n_particulas)
            random = sqrt(rand())
            rand_ang = randn(n_particulas)*2pi
            x[1,:], y[1,:] = condicion_inicial(n_particulas,radio,radio_c) 
            barrera_x, barrera_y = generar_barrera(0,0,25,radio)
            @showprogress "Calculando trayectorias " for i in 2:n_pasos
                


                f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], radio)
                fx, fy = chequear_barrera(x[i-1,:], y[i-1,:],barrera_x,barrera_y,radio)
                #forces1[i-1,:] = f_x
                #forces2[i-1,:] = f_y 
                x[i-1,:] += f_x*dt
                y[i-1,:] += f_y*dt
                x[i-1,:] += fx*dt
                y[i-1,:] += fy*dt


                ruidoDtx = sqrtD * randn(n_particulas)
                
                ruidoDty = sqrtD * randn(n_particulas)
                
                ruidoDr  = sqrtT * randn(n_particulas)
                
                φ[i,:] = φ[i-1,:] .+ Ω*dt +  ruidoDr
                
                x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
                
                y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty

                dx = x[i,:]  .- centro_x
                dy = y[i,:]  .- centro_y
                

                # Verificar la posición de la partícula con respecto al centro del círculo
                #distancia_c = sqrt.((dx).^2 + (dy).^2)

                # Se define un array con el índice de todas las partículas que están fuera de la barrera en el i-ésimo paso temporal
                #fuera = distancia_c .>= (radio_c + radio)  # Tomando en cuenta el radio de la particula

                # Reflejar particulas fuera de la barrera
                #for j in findall(fuera)
                    
                    # Calcular el vector normal al desplazamiento en x e y
                 #   normal_x = dx[j] / distancia_c[j]
                 #   normal_y = dy[j] / distancia_c[j]

                    # Calcular el punto de reflexión tomando en cuenta el radio de la particula
                 #   reflejo_x = centro_x + normal_x * (radio_c + radio)
                 #   reflejo_y = centro_y + normal_y * (radio_c + radio)
                 #   x[i, j] = reflejo_x
                 #   y[i, j] = reflejo_y

                    #Calcular nuevo angulo de orientacion
                 #   angulo = atan(dy[j], dx[j])
                 #   φ[i, j] = angulo + pi / 2
                #end


            end
                

    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",y , ',')
    writedlm("/home/mayron/Datos/barrera_$v/theta_v=00$v.csv",φ , ',')
    writedlm("/home/mayron/Datos/barrera_$v/barrera_x.csv",φ , ',')
    writedlm("/home/mayron/Datos/barrera_$v/barrera_y.csv",φ , ',')
    #writedlm("/home/mayron/Datos/barrera_$v/f_x_v=00$v.csv",forces1 , ',')
    #writedlm("/home/mayron/Datos/barrera_$v/f_y_v=00$v.csv",forces2 , ',')
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
                    # Potencial de interaccion
                    magnitud_fuerza = lj_fuerza(r, 0.5, 2 * radio)
                    # Calculamos la componente x e y de la fuerza    
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
        # Calculate the potential energy
        force = 24 * epsilon * ((sigma^6) / (distancia^8) - 2 * (sigma^12) / (distancia^14))
        return force
end


function condicion_inicial(n_particulas, radio_particula, radio_circulo, max_attempts = 100)
    x_ini = Float64[]  # Array to store x-coordinates
    y_ini = Float64[]  # Array to store y-coordinates



    for _ in 1:n_particulas
        while true
            # Check if the number of attempts exceeds the specified limit

            # Generar condiciones iniciales en coordenadas polares
            angulo = rand() * 2 * π
            r = rand(0:(radio_circulo - radio_particula)) 

            # Convertir a coordenadas cartesianas
            x = r * cos(angulo)
            y = r * sin(angulo)

            # Chequear si la posición de la i-ésima partícula no se solapa con otra
            overlap = false
            for i in 1:length(x_ini)
                if sqrt((x - x_ini[i])^2 + (y - y_ini[i])^2) < 2 * radio_particula
                    overlap = true
                    break
                end
            end

            # Si no hay overlap, entonces guarda la posición en x e y
            if !overlap
                push!(x_ini, x)
                push!(y_ini, y)
                break
            end
        end
    end

    return x_ini, y_ini
end


function quorum_sensing(posicion_x, posicion_y, radio,φ,vc,Ro)
    
    n_particulas = length(posicion_x)
    
    for i in 1:n_particulas
        for j in 1:n_particulas
            if i != j
                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                rij = (dx + dy)/r  

                if rij * (cos(φ[j]) + sin(φ[j])) >= cos(vc) && r <= 4*Ro
                    #Si es que las particulas se encuentran dentro del cono de vision y estan dentro del radio de interaccion
                   #suma = exp^(-rij / Ro)






                    #Aca tendria que calcular la suma de la ec. (3) del paper y luego sumarla a φ[i] ???
                    
                end
            end
        end
    end
    return 
end


function generar_barrera(centro_x, centro_y, radio_c, radio)
    circunferencia = 2 * π * radio_c
    num_particles = round(Int, circunferencia / (2 * radio))

    θ = range(0, stop=2π, length=num_particles+1)[1:end-1]  # Angular positions for particles
    x = centro_x .+ (radio_c - radio) * cos.(θ)
    y = centro_y .+ (radio_c - radio) * sin.(θ)

    return x, y
end


function chequear_barrera(posicion_x, posicion_y, barrera_x, barrera_y, radio)
    n_particulas = length(posicion_x)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
    for i in 1:n_particulas
        for j in 1:length(barrera_x)
            dx = barrera_x[j] - posicion_x[i]
            dy = barrera_y[j] - posicion_y[i]
            r = sqrt(dx^2 + dy^2)  # Distancia entre las particulas y la barrera

            if r <= 2 * radio
                # Potencial de interacción
                magnitud_fuerza = lj_fuerza(r, 0.5, 2 * radio)
                # Calculamos la componente x e y de la fuerza    
                f_x = magnitud_fuerza * dx  
                f_y = magnitud_fuerza * dy  
                # Actualizamos el array x e y de las fuerzas
                fuerza_x[i] += f_x
                fuerza_y[i] += f_y
            end
        end
    end
    return fuerza_x, fuerza_y
end