using Random, DelimitedFiles, ProgressMeter, Distributions, LinearAlgebra, Printf, Dates, Random
folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
L = 500  #diametro del circulo
Dt = 0   #Difusion Traslacional
Dr = 0.01   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)



function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, radio,angulo1::Float64)
   
     # Carpeta donde se guardara los datos de la simulacion
        carpeta = carpeta_simulacion("/home/mayron/Datos")

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, radio, angulo1)
     # Guardar seed
        #Random.seed!(seed)
   
   
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x      = zeros(n_pasos,n_particulas)
        y      = similar(x)
        φ      = similar(x)
        #m      = Vector{Float64}[] #magnetizacion
        x_data = Vector{Float64}[]
        y_data = Vector{Float64}[]
        φ[1,:] = rand(0:2pi,n_particulas).*randn(n_particulas)
        x[1,:] , y[1,:] = condicion_inicial(n_particulas,radio,48)
        push!(x_data, x[1, :])
        push!(y_data, y[1, :])
        #m[1] = abs.((sum( ( cos.( φ[1,:] ) ) + ( sin.( φ[1,:] ) ) ))/n_particulas)
        barrera_x, barrera_y = generar_barrera(0,0,50,radio)
        @showprogress "Calculando trayectorias " for i in 2:n_pasos
            

            
            f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], radio)
            quorum, Nc = quorum_sensing(x[i-1,:],y[i-1,:],n_particulas,φ[i-1,:],angulo1)
            τ = torque_barrera(x[i-1,:], y[i-1,:],barrera_x,barrera_y,φ[i-1,:], radio)
            x[i-1,:] += f_x * dt
            y[i-1,:] += f_y * dt




            ruidoDtx = sqrtD * randn(n_particulas)
            
            ruidoDty = sqrtD * randn(n_particulas)
            
            ruidoDr  = sqrtT * randn(n_particulas)
            
            φ[i,:] = φ[i-1,:] + 5*(quorum./Nc)*dt   +  τ./100   +  ruidoDr
            
            x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
            
            y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty

            #m[i] = abs.(sum( ( cos.( φ[i,:] ) ) + ( sin.( φ[i,:] ) ) )/(n_particulas))


            #x[i,:], y[i,:] = periodic_bc(x[i,:],y[i,:],n_particulas,L)

            
                if i % 100 == 0
                    # Save data for every 100th step
                    push!(x_data, x[i, :])
                    push!(y_data, y[i, :])
                end

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
      
    return x_data, y_data
end


function correccion_lj(posicion_x, posicion_y, radio)
    n_particulas = length(posicion_x)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
    Threads.@threads for i in 1:n_particulas
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


function quorum_sensing(posicion_x, posicion_y, n_particulas, φ, angulo1, Ro = 3)
    
    quorum = zeros(n_particulas)
    Nc = ones(n_particulas)
    
    Threads.@threads for i in 1:n_particulas
        for j in 1:n_particulas
            if i != j

                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                rij = [dx, dy] / r  


                #angulo = atan(rij[2],rij[1])

                if (dot(rij, [cos(φ[i]), sin(φ[i])]) >= cos(angulo1)) && (r <= 4 * Ro)
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
    Threads.@threads for i in 1:n_particulas
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



function generar_barrera(centro_x, centro_y, radio_c, radio=0.5)
    circunferencia = 2 * π * radio_c
    num_particles = round(Int, circunferencia / (2 * radio))

    θ = range(0, stop=2π, length=num_particles+1)[1:end-1]  # Angular positions for particles
    x = centro_x .+ (radio_c - radio) * cos.(θ)
    y = centro_y .+ (radio_c - radio) * sin.(θ)

    return x, y
end


function torque_barrera(posicion_x, posicion_y, barrera_x, barrera_y, φ, radio)
    n_particulas = length(posicion_x)
    torque = zeros(n_particulas)
    
    for i in 1:n_particulas
        Threads.@threads for j in 1:length(barrera_x)
            dx = barrera_x[j] - posicion_x[i]
            dy = barrera_y[j] - posicion_y[i]
            r = sqrt(dx^2 + dy^2)  # Distancia entre las particulas y la barrera

            if r <= 2 * radio * 2^(1/6)
                # Calculamos "n_wall"
                N_wall = [barrera_x[j], barrera_y[j], 0] / norm([barrera_x[j], barrera_y[j]])
                t_wall = cross(N_wall, [0, 0, 1])

                # Calculo del torque
                torque[i] = -(dot([cos(φ[i]), sin(φ[i]), 0], N_wall)) * (dot([cos(φ[i]), sin(φ[i]), 0], t_wall))
            end
        end
    end
    
    return torque
end

function carpeta_simulacion(base_dir)
    #Obtener la fecha de hoy
    fecha = Dates.today()

    #Contador de las simulaciones realizadas en el dia
    contador_sim = 1

    #Formatear la fecha para darle nombre a la carpeta
    nombre_carpeta = Dates.format(fecha, "yyyy-dd-mm")

    #Si ya existe una carpeta en el dia, le añade un numero como sufijo para no sobreescribir
    while isdir(joinpath(base_dir, "$nombre_carpeta-$contador_sim"))
        contador_sim += 1
    end

    #Se crea la carpeta
    nombre_carpeta = joinpath(base_dir, "$nombre_carpeta-$contador_sim")
    mkdir(nombre_carpeta)

    return nombre_carpeta
end

function generar_log(folder_path, v, n_pasos, n_particulas, radio, angulo1)
    # Formatear la velocidad para incluirlo como string
    v_str = @sprintf("%.2f", v)

    # Escribir los parámetros en el archivo log
    log_filename = joinpath(folder_path, "log.txt")
    open(log_filename, "w") do file
        println(file, "ESta simulación se realizó el: $(Dates.now())")
        println(file, "Parámetros de la simulación:")
        println(file, "v: $v")
        println(file, "n_pasos: $n_pasos")
        println(file, "n_particulas: $n_particulas")
        println(file, "radio: $radio")
        println(file, "angulo1: $angulo1")
        #println(file, "Seed: $seed")
        println(file, "--------------------------------------")
    end
end


