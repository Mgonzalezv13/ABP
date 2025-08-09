using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-4  #Paso temporal




function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, L::Int64, angulo1::Float64,α,Dr,radio_p = 1)
   
     # Carpeta donde se guardara los datos de la simulacion
     carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α)
        η       = packing_fraction(n_particulas,L)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, L, angulo1,η,α,Dr,dt)
     # Guardar seed
        #Random.seed!(seed)
   
        sqrtT = sqrt(2*Dr*dt)
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x_data   = Vector{Float64}[]
        y_data   = Vector{Float64}[]
        vx_data  = Vector{Float64}[]
        vy_data  = Vector{Float64}[]
        φ_data   = Vector{Float64}[]
        φ_old = rand(0:2pi,n_particulas)
        radio_p = rand(radio_p:0.4:1.4,n_particulas)
        x_old , y_old = ini_con_pbc(n_particulas,L,radio_p)
        vf,vc = vecinos_pbc(x_old,y_old,α,L)
        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y = correccion_lj(x_old, y_old, vf,radio_p,n_particulas)
            quorum, Nc = quorum_sensing(x_old, y_old, n_particulas,φ_old, angulo1,vc,α)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)

            vx = v*cos.(φ_old)  + f_x
            
            vy = v*sin.(φ_old)  + f_y

            φ = φ_old + 5*(quorum./Nc)*dt + ruidoDr 
            
            x = x_old + vx*dt
            
            y = y_old + vy*dt

            #reflexion
            x, y = periodic_bc(x,y,n_particulas,L)  
            
                if i % 100 == 0
                    # Guardar datos cada 100 pasos de tiempo
                    push!(x_data, x)
                    push!(y_data, y)
                    push!(φ_data, φ)
                    push!(vx_data,vx)
                    push!(vy_data,vy)
                end



               
                
                if i % 200 == 0

                    #actualizar la lista de vecinos cada ciertos pasos
                    vf,vc = vecinos_pbc(x,y,α,L)

                end
                
        

                


                φ_old = φ
                x_old = x
                y_old = y 

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
        writedlm(joinpath(carpeta, "vx=$(round(v, digits=2)).csv"), vx_data, ',')
        writedlm(joinpath(carpeta, "vy=$(round(v, digits=2)).csv"), vy_data, ',')
        writedlm(joinpath(carpeta, "radio=$(round(v, digits=2)).csv"), radio_p, ',')
        
    return
end


function correccion_lj(posicion_x, posicion_y,vecinos, radio,n_particulas)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
     for i in 1:n_particulas
        for j in vecinos[i]
                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]

                dx -= L * round(dx / L)
                dy -= L * round(dy / L)
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                sigma = (radio[j] + radio[i])/2
                # Potencial de interaccion
                if r <= 2^(1/6)*2*sigma
                    magnitud_fuerza = lj_fuerza(r, 0.5, 2 * sigma)
                    # Calculamos la componente x e y de la fuerza    
                    f_x = magnitud_fuerza * dx 
                    f_y = magnitud_fuerza * dy 
                    # Updateamos el array x e y de las fuerzas
                    fuerza_x[i] += f_x
                    fuerza_y[i] += f_y
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


function condicion_inicial(radio_circulo::Int64, radio_particula::Vector{Float64}, n_particulas::Int64; max_attempts::Int = 1000)
    x_ini = Float64[]  # Array to store x-coordinates
    y_ini = Float64[]  # Array to store y-coordinates

    for j in 1:n_particulas
        intento = 0
        colocado = false

        while intento < max_attempts
            intento += 1

            # Generar condiciones iniciales en coordenadas polares
            angulo = rand() * 2π
            r = rand() * (radio_circulo - radio_particula[j])  # Evita colocarla fuera del círculo

            # Convertir a coordenadas cartesianas
            x = r * cos(angulo)
            y = r * sin(angulo)

            # Chequear si la partícula no se solapa con otra
            overlap = false
            for i in 1:length(x_ini)
                dx = x - x_ini[i]
                dy = y - y_ini[i]
                distancia = sqrt(dx^2 + dy^2)
                if distancia < (radio_particula[i] + radio_particula[j])
                    overlap = true
                    break
                end
            end

            # Si no hay overlap, guarda la posición
            if !overlap
                push!(x_ini, x)
                push!(y_ini, y)
                colocado = true
                break
            end
        end

        if !colocado
            error("No se pudo colocar la partícula $j después de $max_attempts intentos.")
        end
    end

    return x_ini, y_ini
end



function quorum_sensing(posicion_x, posicion_y, n_particulas, φ, angulo1,vecinos,α,Ro=3)
    
    quorum = zeros(n_particulas)
    Nc = ones(n_particulas)
    
     for i in 1:n_particulas
        for j in vecinos[i]
            if i != j
                
                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]

                dx -= L * round(dx / L)
                dy -= L * round(dy / L)

                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
                rij = [dx, dy] / r      


                #angulo = atan(rij[2],rij[1])

                if (dot(rij, [cos(φ[i]), sin(φ[i])]) >= cos(angulo1)) && (r <= α * Ro)
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






function carpeta_simulacion(base_dir, angulo1, n_particulas, Dr, α)
    # Asegurarse que el directorio base exista
    mkpath(base_dir)

    # Representación bonita del ángulo (usar isapprox para evitar problemas float)
    angle_str = if isapprox(angulo1, π)
        "π"
    elseif isapprox(angulo1, π/2)
        "π_2"
    elseif isapprox(angulo1, π/3)
        "π_3"
    elseif isapprox(angulo1, π/4)
        "π_4"
    elseif isapprox(angulo1, π/10)
        "π_10"
    else
        string(angulo1)
    end

    # Nombre base con parámetros
    nombre_base = "PBC_Poli_N=$(n_particulas)-Dr=$(Dr)-θ=$(angle_str)-Ro=$(3*α)"

    # Construir nombre único: si ya existe, añadir _2, _3, ...
    contador = 1
    candidato = joinpath(base_dir, nombre_base)
    while isdir(candidato) || isfile(candidato)
        contador += 1
        candidato = joinpath(base_dir, "$(nombre_base)_$contador")
    end

    # Crear la carpeta y devolver la ruta completa
    mkdir(candidato)
    return candidato
end

function generar_log(folder_path, v, n_pasos, n_particulas, R, angulo1, η, α, Dr, dt)
    log_filename = joinpath(folder_path, "log.txt")
    open(log_filename, "w") do file
        println(file, "Parámetros de la simulación:")
        println(file, "Ángulo: $angulo1")
        println(file, "N partículas: $n_particulas")
        println(file, "Dr: $Dr")
        println(file, "dt: $dt")
        println(file, "v: $v")
        println(file, "Iteraciones: $n_pasos")
        println(file, "Radio Confinamiento: $R")
        println(file, "Empaquetamiento: $η")
        println(file, "Tamaño cono: $(3*α)")
    end
end




function reflective_bc(x_old, y_old, x, y, R)
    N = length(x_old)
    r = @. sqrt(x^2 + y^2)          # Distance from origin
    outside_mask = r .>= R           # Mask for particles outside the circle
    outside_indices = findall(outside_mask)

    # Preallocate intersection points
    p_x = zeros(N)
    p_y = zeros(N)

    # Compute intersection points for particles outside
    for i in outside_indices
        m = (y[i] - y_old[i]) / (x[i] - x_old[i])  # Slope of trajectory line
        h = y_old[i] - m * x_old[i]                # y-intercept
        mm = m^2
        delta = sqrt(R^2 * (1 + mm) - h^2)          # Discriminant
        raiz1_x = (-m * h + delta) / (1 + mm)       # First intersection (x)
        raiz1_y = m * raiz1_x + h                   # First intersection (y)
        raiz2_x = (-m * h - delta) / (1 + mm)       # Second intersection (x)
        raiz2_y = m * raiz2_x + h                   # Second intersection (y)

        # Choose the correct intersection (between old and new positions)
        if (x_old[i] < raiz1_x < x[i]) || (y_old[i] < raiz1_y < y[i]) ||
           (x[i] < raiz1_x < x_old[i]) || (y[i] < raiz1_y < y_old[i])
            p_x[i], p_y[i] = raiz1_x, raiz1_y
        else
            p_x[i], p_y[i] = raiz2_x, raiz2_y
        end
    end

    # Compute normal vectors and reflections
    n_x = @. (-1 / R) * p_x * outside_mask  # Normal vector (inward)
    n_y = @. (-1 / R) * p_y * outside_mask
    factor = @. (x - p_x) * n_x + (y - p_y) * n_y  # Dot product
    x_reflected = @. (x - 2 * n_x * factor) * outside_mask
    y_reflected = @. (y - 2 * n_y * factor) * outside_mask

    # Combine results: keep particles inside, reflect those outside
    x_final = @. x * (r < R) + x_reflected
    y_final = @. y * (r < R) + y_reflected

    return x_final, y_final
end

function packing_fraction(N, L, r=1)
    η = (π * N * r^2) / L^2
    return η
end




function vecinos_pbc(x, y, rq,L)
    N = length(x)
    inv_L = 1.0 / L

    veci_fuerza = Vector{Int}[]
    veci_cono = Vector{Int}[]

    r_fuerza = 2*2^(1/6) +1.0
    r_quorum = rq*3 +1.0
   
    x1 = x .- x'
    y1 = y .- y'

    dx = x1 .- L .* round.(x1 .* inv_L)
    dy = y1 .- L .* round.(y1 .* inv_L)
    
    r_min = sqrt.(dx.^2 + dy.^2)

    
    for i in 1:size(r_min,2) 

        v_fuerza = findall(0 .< r_min[i,:] .<= r_fuerza)

        v_quorum = findall(0 .< r_min[i,:] .<= r_quorum)

        push!(veci_fuerza, v_fuerza)

        push!(veci_cono, v_quorum)

    end
    
    return veci_fuerza, veci_cono
end




function ini_con_pbc(N, L, radio_particula::Vector{Float64})

    x_ini = Float64[]
    y_ini = Float64[]

    for n in 1:N
        attempts = 0
        while true
            attempts += 1

            # Random position in [-L/2, L/2]
            x = rand() * L - L/2
            y = rand() * L - L/2

            # Check overlaps with PBC
            overlap = false
            for i in 1:length(x_ini)
                dx = abs(x - x_ini[i])
                dy = abs(y - y_ini[i])

                # Apply periodic boundary conditions
                dx = min(dx, L - dx)
                dy = min(dy, L - dy)

                # Minimum allowed distance = sum of radii
                min_dist = radio_particula[n] + radio_particula[i]

                if sqrt(dx^2 + dy^2) < min_dist
                    overlap = true
                    break
                end
            end

            if !overlap
                push!(x_ini, x)
                push!(y_ini, y)
                break
            end
        end
    end

    return x_ini, y_ini
end
