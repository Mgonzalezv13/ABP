using Random, DelimitedFiles, ProgressMeter, Distributions, LinearAlgebra

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 100  #diametro del circulo
Dt = 0   #Difusion Traslacional
Dr = 0.01   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-4  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)


function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, radio)
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x   = zeros(n_pasos,n_particulas)
        y   = similar(x)
        φ   = similar(x)
        quorum1 = zeros(n_pasos,n_particulas)
         φ[1,:] = rand(0:2pi,n_particulas)
         x[1,:] = rand(-30:30,n_particulas)
         y[1,:] = rand(-30:30,n_particulas)
        @showprogress "Calculando trayectorias " for i in 2:n_pasos
            


            f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], radio)
            quorum, Nc = quorum_sensing(x[i-1,:],y[i-1,:],n_particulas,φ[i-1,:])
            quorum1[i-1,:] = quorum
            x[i-1,:] += f_x*dt
            y[i-1,:] += f_y*dt
            #forces1[i-1,:] = f_x
            #forces2[i-1,:] = f_y 




            ruidoDtx = sqrtD * randn(n_particulas)
            
            ruidoDty = sqrtD * randn(n_particulas)
            
            ruidoDr  = sqrtT * randn(n_particulas)
            
            φ[i,:] = 5*(quorum./Nc)*dt + φ[i-1,:] +  ruidoDr
            
            x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
            
            y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty


            x[i,:], y[i,:] = periodic_bc(x[i,:],y[i,:],n_particulas,L)


            end  

    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",y , ',')
    writedlm("/home/mayron/Datos/barrera_$v/theta_v=00$v.csv",φ , ',')
   # writedlm("/home/mayron/Datos/barrera_$v/fuerza_x.csv",forces1 , ',')
   # writedlm("/home/mayron/Datos/barrera_$v/fuerza_y.csv",forces2 , ',')
    writedlm("/home/mayron/Datos/barrera_$v/quorum.csv",quorum1 , ',')
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


function quorum_sensing(posicion_x, posicion_y, n_particulas, φ, Ro = 3)
    
    quorum = zeros(n_particulas)
    Nc = ones(n_particulas)
    
    for i in 1:n_particulas
        for j in 1:n_particulas
            if i != j

                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                rij = [dx, dy] / r  


                #angulo = atan(rij[2],rij[1])

                if (dot(rij, [cos(φ[i]), sin(φ[i])]) >= cos(pi/3)) && (r <= 4 * Ro)
                    # si las particulas estan dentro del cono de vision y a la distancia indicada


                    angulo = atan(rij[2],rij[1])

                    quorum[i] += exp(-r/ Ro)*sin(angulo - φ[i])
                    Nc[i] += exp(-r/ Ro)
                end
            end
        end
    end 
    return quorum, Nc
end


function periodic_bc(posicion_x, posicion_y, n_particulas, L)
    for i in 1:n_particulas
        if posicion_x[i] > L/2
            posicion_x[i] -= L
        end

        if posicion_x[i] < -L/2
            posicion_x[i] += L
        end

        if posicion_y[i] > L/2
            posicion_y[i] -= L
        end

        if posicion_y[i] < -L/2
            posicion_y[i] += L
        end
    end

    return posicion_x, posicion_y
end
