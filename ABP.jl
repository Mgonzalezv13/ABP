using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-4  #Paso temporal




function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, R::Float64, angulo1::Float64,α,Dr,radio_p = 1)
   
     # Carpeta donde se guardara los datos de la simulacion
        carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α)
        η       = packing_fraction(n_particulas,R)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, R, angulo1,η,α,Dr,dt)
    
        sqrtT = sqrt(2*Dr*dt) #esta cantidad se mantiene fija

        x_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        y_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        vx_data  = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        vy_data  = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        φ_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)

        φ_old = rand(0:2pi,n_particulas)
        x_old , y_old = condicion_inicial(n_particulas,R)
        vx_old , vy_old = zeros(n_particulas), zeros(n_particulas) 
        vf,vc = vecinos(x_old,y_old,α)
        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y = correccion_lj(x_old, y_old, vf,radio_p,n_particulas)
            quorum, Nc = quorum_sensing(x_old, y_old, n_particulas,φ_old, angulo1,vc,α)
            
            #tomando γ = m = 1

            γ = 1.0     
            m = 1.0     
            
            ruidoDr  = sqrtT * randn(n_particulas)

            vx = vx_old + dt/m * ( -γ*vx_old + γ*v*cos.(φ_old)  + f_x )

            vy = vy_old + dt/m * ( -γ*vy_old + γ*v*sin.(φ_old)  + f_y )
            
            φ  = φ_old + ruidoDr + 5*(quorum./Nc)*dt  
            
            x  = x_old + vx*dt
            
            y  = y_old + vy*dt

            x, y = reflective_bc(x_old,y_old,x,y,R) 
            
            
                if i % 100 == 0
                    # Guardar datos cada 100 pasos de tiempo
                    x_data[Int64(i/100),:]  .= x 
                    y_data[Int64(i/100),:]  .= y
                    φ_data[Int64(i/100),:]  .= φ
                    vx_data[Int64(i/100),:] .= vx
                    vy_data[Int64(i/100),:] .= vy
                end



               
                
                if i % 200 == 0

                    #actualizar la lista de vecinos cada ciertos pasos
                    vf,vc = vecinos(x, y, α)

                end
                
        
 


                φ_old = φ
                x_old = x
                y_old = y 
                vx_old = vx
                vy_old = vy

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
        writedlm(joinpath(carpeta, "vx=$(round(v, digits=2)).csv"), vx_data, ',')
        writedlm(joinpath(carpeta, "vy=$(round(v, digits=2)).csv"), vy_data, ',')

    return
end

function correccion_lj(posicion_x, posicion_y,vecinos, radio,n_particulas)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)
    
     for i in 1:n_particulas
        for j in vecinos[i]

                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
                r = sqrt(dx^2 + dy^2)  # Distancia entre la i-esima y j-esima particula
               
                if r <= 2*radio
                    # Potencial de interaccion
                    magnitud_fuerza = soft_fuerza(r, 100, 2 * radio)
                    # Calculamos la componente x e y de la fuerza    
                    f_x = -magnitud_fuerza * (dx/r) 
                    f_y = -magnitud_fuerza * (dy/r) 
                    # Updateamos el array x e y de las fuerzas
                    fuerza_x[i] += f_x
                    fuerza_y[i] += f_y
                end
        end
    end
    return fuerza_x, fuerza_y
end

function soft_fuerza(distancia, epsilon, sigma)
    return epsilon * (1 - (distancia /(sigma)))^(3/2)
end


function condicion_inicial(n_particulas,radio_circulo,radio_particula = 1, max_attempts = 100)
    x_ini = Float64[]  # Array to store x-coordinates
    y_ini = Float64[]  # Array to store y-coordinates



    for _ in 1:n_particulas
        while true

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


function quorum_sensing(posicion_x, posicion_y, n_particulas, φ, angulo1,vecinos,α,Ro=3)
    
    quorum = zeros(n_particulas)
    Nc = ones(n_particulas)
    
     for i in 1:n_particulas
        for j in vecinos[i]
            if i != j

                dx = posicion_x[j] - posicion_x[i]
                dy = posicion_y[j] - posicion_y[i]
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

function periodic_bc(posicion_x, posicion_y, L)
    posicion_x = mod.(posicion_x .+ L/2, L) .- L/2
    posicion_y = mod.(posicion_y .+ L/2, L) .- L/2
    return posicion_x, posicion_y
end




function carpeta_simulacion(base_dir, angulo1, n_particulas, Dr,α)
    # Convert angle to nice π format if it's a multiple of π
    angle_str = if angulo1 == π
        "π"
    elseif angulo1 == π/2
        "π_2"
    elseif angulo1 == π/3
        "π_3"
    elseif angulo1 == π/4
        "π_4"
    elseif angulo1 == π/10
        "π_10"
    else
        "$angulo1"
    end

    # Create base folder name with parameters
    param_str = "Barrera_N=$(n_particulas)-Dr=$(Dr)-θ=$(angle_str)-Ro=$(3*α)"
    
    # Check if folder exists
    contador_sim = 1
    nombre_carpeta = joinpath(base_dir, param_str)
    
    # If exists, add repeticion counter
    if isdir(nombre_carpeta)
        while isdir(joinpath(base_dir, "$param_str-repeticion_$contador_sim"))
            contador_sim += 1
        end
        nombre_carpeta = joinpath(base_dir, "$param_str-repeticion_$contador_sim")
    end
    
    mkdir(nombre_carpeta)
    return nombre_carpeta
end

function generar_log(folder_path, v, n_pasos, n_particulas, R, angulo1, η, α, Dr, dt)
    # Write parameters to log file
    log_filename = joinpath(folder_path, "log.txt")
    open(log_filename, "w") do file
        println(file, "Parámetros de la simulación:")
        println(file, "Ángulo: $angulo1")
        println(file, "N partículas: $n_particulas")
        println(file, "Dr: $Dr")
        println(file, "dt: $dt")
        println(file, "v: $v")
        println(file, "Iteraciones: $n_pasos")
        println(file, "Radio Barrera: $R")
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

function packing_fraction(N, R, r=1)
    η = (N * r^2) / R^2
    return η
end




function vecinos(x,y, rq = 2)
    
    veci_fuerza = Vector{Int}[]
    veci_cono = Vector{Int}[]
    pos = hcat(x,y)


    r = pairwise(Euclidean(), pos, dims=1)

    r_fuerza = 2*2^(1/6) +1.0
    r_quorum = rq*3 +1.0

    for i in 1:size(r,2) 

        v_fuerza = findall(0 .< r[i,:] .<= r_fuerza)

        v_quorum = findall(0 .< r[i,:] .<= r_quorum)

        push!(veci_fuerza, v_fuerza)

        push!(veci_cono, v_quorum)

    end


    return veci_fuerza, veci_cono
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

function ini_con_pbc(N, L, radio_particula=1)
    x_ini = Float64[]
    y_ini = Float64[]
    min_dist = 2 * radio_particula

    for _ in 1:N
        attempts = 0
        while true
            attempts += 1

            # Continuous random positions instead of grid
            x = rand() * L - L/2  # Random position in [-L/2, L/2]
            y = rand() * L - L/2

            # Check overlaps with PBC
            overlap = false
            for i in 1:length(x_ini)
                dx = abs(x - x_ini[i])
                dy = abs(y - y_ini[i])
                # Apply periodic boundary conditions
                dx = min(dx, L - dx)
                dy = min(dy, L - dy)
                
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

