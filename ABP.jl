using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-3  #Paso temporal




function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, radio::Int64, angulo1::Float64,α,Dr,radio_p = 1, radio_b = 0.1)
   
     # Carpeta donde se guardara los datos de la simulacion
     carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α)
        η       = packing_fraction(n_particulas,radio)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, radio, angulo1,η,α,Dr,dt)
     # Guardar seed
        #Random.seed!(seed)
   
        sqrtT = sqrt(2*Dr*dt)
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        y_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        vx_data  = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        vy_data  = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        φ_data   = Matrix{Float64}(undef, Int64(n_pasos/100), n_particulas)
        φ_old = rand(0:2pi,n_particulas)
        x_old , y_old = condicion_inicial(n_particulas,radio)
        b_x, b_y = generar_barrera(0,0,radio,radio_b)
        vf,vc = vecinos(x_old,y_old,α)
        vb    = vecinos_barrera(x_old,y_old,b_x,b_y)
        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y    = correccion_lj(x_old, y_old, vf,radio_p,n_particulas)
            f_bx, f_by  = correccion_barrera(x_old, y_old, b_x, b_y,vb, radio_p,radio_b, n_particulas)
            quorum, Nc  = quorum_sensing(x_old, y_old, n_particulas,φ_old, angulo1,vc,α)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)

            vx = v*cos.(φ_old)  + f_x + f_bx
            
            vy = v*sin.(φ_old)  + f_y + f_by

            φ = φ_old + 5*(quorum./Nc)*dt + ruidoDr 
            
            x = x_old + vx*dt
            
            y = y_old + vy*dt

            
            
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
                    vf,vc = vecinos(x,y,α)
                    vb    = vecinos_barrera(x,y,b_x,b_y)

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
        
    return
end


function correccion_barrera(pos_x, pos_y, barrera_x, barrera_y,vecinos_b, radio_p,radio_bar, n_particulas)

    fx = zeros(n_particulas)
    fy = zeros(n_particulas)

    for i in 1:n_particulas
        for j in vecinos_b[i]

        dx = barrera_x[j] - pos_x[i]
        dy = barrera_y[j] - pos_y[i]
        r = sqrt(dx^2 + dy^2)

        sigma = (radio_bar + radio_p)/2

            if r <= 2 * sigma
            F = soft_fuerza(r, 125, 2 * sigma)

            fx[i] -= F * (dx / r)
            fy[i] -= F * (dy / r)
            end
        end
    end
        return fx, fy
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
                    magnitud_fuerza = soft_fuerza(r, 75, 2 * radio)
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

    
    param_str = "Barrera_N=$(n_particulas)-θ=$(angle_str)-Dr=$(Dr)-Ro=$(3*α)"
    nombre_carpeta = joinpath(base_dir, param_str)

    
    contador_sim = 1
    nombre_carpeta = joinpath(base_dir, param_str)
    
    
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
        println(file, "Radio Confinamiento: $R")
        println(file, "Empaquetamiento: $η")
        println(file, "Tamaño cono: $(3*α)")
    end
end


function packing_fraction(N, R, r=1)
    η = N * r^2 / R^2
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


function generar_barrera(centro_x, centro_y, radio_c, radio=0.5)
    circunferencia = 2 * π * radio_c
    num_particles = round(Int, circunferencia / (2 * radio))

    θ = range(0, stop=2π, length=num_particles+1)[1:end-1]  # Angular positions for particles
    x = centro_x .+ (radio_c - radio) * cos.(θ)
    y = centro_y .+ (radio_c - radio) * sin.(θ)

    return x, y
end



function vecinos_barrera(x_spp, y_spp, x_bar, y_bar, cutoff=3.0)
    n_spp = length(x_spp)
    n_bar = length(x_bar)

    vec = Vector{Vector{Int}}(undef, n_spp)

    for i in 1:n_spp
        list = Int[]
        for j in 1:n_bar
            dx = x_bar[j] - x_spp[i]
            dy = y_bar[j] - y_spp[i]
            if dx^2 + dy^2 <= cutoff^2
                push!(list, j)
            end
        end
        vec[i] = list
    end

    return vec
end