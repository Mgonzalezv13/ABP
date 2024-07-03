using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics


Dt = 0   #Difusion Traslacional
Dr = 8e-2   #Difusion Rotacional
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
        x       = Matrix{Float64}(undef,n_pasos,n_particulas)
        y       = Matrix{Float64}(undef,n_pasos,n_particulas)
        φ       = Matrix{Float64}(undef,n_pasos,n_particulas)
        p       = Float64[]
        rg      = Float64[]
        x_data  = Vector{Float64}[]
        y_data  = Vector{Float64}[]
        φ_data  = Vector{Float64}[]
        φ[1,:] = rand(0:2pi,n_particulas)
        x[1,:] , y[1,:] = condicion_inicial(n_particulas,radio,48)
        barrera_x, barrera_y = generar_barrera(0,0,50,radio)
        vecinos = Verlet_vecinos(n_particulas,x[1,:],y[1,:],13.)
        v_barr = Verlet_vecinos(n_particulas,x[1,:],y[1,:],barrera_x,barrera_y,5.)
        vecinos_q = Verlet_vecinos(n_particulas,x[1,:],y[1,:],16.)
         for i in 2:n_pasos

          
            f_x, f_y = correccion_lj(x[i-1,:], y[i-1,:], vecinos,radio,n_particulas)
            fb_x, fb_y = chequear_barrera(x[i-1,:], y[i-1,:],barrera_x,barrera_y, radio, v_barr, n_particulas)
            quorum, Nc = quorum_sensing(x[i-1,:],y[i-1,:],n_particulas,φ[i-1,:], angulo1)
            τ = torque_barrera(x[i-1,:], y[i-1,:],barrera_x,barrera_y,φ[i-1,:],  radio, v_barr,n_particulas)
            x[i-1,:] += f_x * dt
            y[i-1,:] += f_y * dt
            x[i-1,:] += fb_x*dt
            y[i-1,:] += fb_y*dt



            ruidoDtx = sqrtD * randn(n_particulas)
            
            ruidoDty = sqrtD * randn(n_particulas)
            
            ruidoDr  = sqrtT * randn(n_particulas)
            
            φ[i,:] = φ[i-1,:] + 5*(quorum./Nc)*dt   +  τ./10   +  ruidoDr
            
            x[i,:] = x[i-1,:] + v*cos.(φ[i-1,:])*dt + ruidoDtx
            
            y[i,:] = y[i-1,:] + v*sin.(φ[i-1,:])*dt +  ruidoDty
            
            
                if i % 100 == 0
                    # Guardar datos cada 100 pasos de tiempo
                    push!(x_data, x[i, :])
                    push!(y_data, y[i, :])
                    push!(φ_data, φ[i, :])
                    # Calcular el radio de giro y la magnetizacion cada 100 pasos
                    push!(p,abs.(sum( ( cos.( φ[i,:] ) ) + ( sin.( φ[i,:] ) ) )/(n_particulas)))
                    push!(rg,rg_dt(x[i,:],y[i,:],10,n_particulas))

                end




                if i % 200 == 0
                    #actualizar la lista de vecinos cada 200 pasos
                    vecinos = Verlet_vecinos(n_particulas,x[i,:],y[i,:],13.)
                    v_barr = Verlet_vecinos(n_particulas,x[i,:],y[i,:],barrera_x,barrera_y,5.)
                    vecinos_q = Verlet_vecinos(n_particulas,x[i,:],y[i,:],16.)
        
                end

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
        writedlm(joinpath(carpeta, "rg_N=$n_particulas.csv"), rg, ',')
        writedlm(joinpath(carpeta, "P_N=$n_particulas.csv"), p, ',')

        
      
    return p,rg
end


function correccion_lj(posicion_x, posicion_y,vecinos, radio,n_particulas)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
     for i in 1:n_particulas
        for j in vecinos[i]
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


function quorum_sensing(posicion_x, posicion_y, n_particulas, φ, angulo1,Ro = 3)
    
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



function generar_barrera(centro_x, centro_y, radio_c, radio=0.5)
    circunferencia = 2 * π * radio_c
    num_particles = round(Int, circunferencia / (2 * radio))

    θ = range(0, stop=2π, length=num_particles+1)[1:end-1]  # Angular positions for particles
    x = centro_x .+ (radio_c - radio) * cos.(θ)
    y = centro_y .+ (radio_c - radio) * sin.(θ)

    return x, y
end

function chequear_barrera(posicion_x, posicion_y, barrera_x, barrera_y, radio, vecinos, n_particulas)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    for i in 1:n_particulas
        for j in vecinos[i]
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

function torque_barrera(posicion_x, posicion_y, barrera_x, barrera_y, φ, radio, vecinos,n_particulas)
    torque = zeros(n_particulas)
    
    for i in 1:n_particulas
         for j in vecinos[i]
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



function rg_dt(x::Vector{Float64}, y::Vector{Float64}, k_clusters,n_particulas)
    
    rg = Vector{Float64}(undef,n_particulas)

   
    pos = hcat(x, y)

    
    cluster_a = kmeans(pos', k_clusters)

   
    ind_part = cluster_a.assignments

    for i in 1:n_particulas
        # Veo en que cluster esta mi particula e identifico el indice de ese cluster
        ind_cluster = ind_part[i]

        # distancia de la i-esima particula al centro del cluster 
        distancia= sqrt(sum((pos[i,:] .- cluster_a.centers[:,ind_cluster]).^2))

        # Calculate radius of gyration for the ith particle
        rg[i] = distancia
    end

    r_g = mean(rg)



    return r_g
end

function Verlet_vecinos(n_particulas::Int64, posicion_x::Vector{Float64}, posicion_y::Vector{Float64}, r_interaccion::Float64, delta_r=1)
    # Define un array donde se guardan otro array con los indices de los vecinos de la particula i
    veci_verlet = Vector{Int64}[]
    # Defino un radio efectivo para chequear los vecinos
    r_efectivo = r_interaccion + delta_r

    for i in 1:n_particulas
        vecinos_i = Int[]

        for j in 1:n_particulas
            if i != j

                dx  = posicion_x[j] - posicion_x[i]
                dy  = posicion_y[j] - posicion_y[i]
                r   = sqrt(dx^2 + dy^2)

                if r < r_efectivo
                    push!(vecinos_i, j)
                end
            end    
        end

        push!(veci_verlet, vecinos_i)
    end

    return veci_verlet
end

function Verlet_vecinos(n_particulas::Int64, posicion_x::Vector{Float64}, posicion_y::Vector{Float64}, barrera_x::Vector{Float64}, barrera_y::Vector{Float64}, r_interaccion::Float64, delta_r=1)
    
    veci_verlet = Vector{Int64}[]
 
    r_efectivo = r_interaccion + delta_r

    for i in 1:n_particulas
        vecinos_i = Int[]

        for j in 1:length(barrera_x)
            dx = barrera_x[j] - posicion_x[i]
            dy = barrera_y[j] - posicion_y[i]
            r = sqrt(dx^2 + dy^2)

            if r < r_efectivo
                push!(vecinos_i, j)
            end
        end    

        push!(veci_verlet, vecinos_i)
    end

    return veci_verlet
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


