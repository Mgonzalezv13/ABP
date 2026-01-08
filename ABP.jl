using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-4  #Paso temporal

intervalo = 2000


function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, L::Float64, angulo1::Float64,α::Float64,Dr::Float64,radio_p = 1)
   
     # Carpeta donde se guardara los datos de la simulacion
        carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α)
        η       = packing_fraction(n_particulas,L)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, L, angulo1,η,α,Dr,dt)
    
        sqrtT = sqrt(2*Dr*dt) #esta cantidad se mantiene fija

        x_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        y_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vx_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vy_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        φ_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)

        #x_old , y_old, φ_old = ini_pbc(n_particulas,L,α)
        vx, vy = zeros(n_particulas), zeros(n_particulas)
        x, y, φ = ini_pbc(n_particulas,L,α)
        vf, vc, dx, dy, r2 = vecinos_pbc(x,y,L,α)
        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y = correccion_lj(vf,dx,dy,r2,radio_p,n_particulas)
            quorum, Nc = quorum_sensing(φ,vc,dx,dy,r2,angulo1,α,3*radio_p)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)

            @. vx = v*cos(φ)*dt + f_x*dt
            @. vy = v*sin(φ)*dt + f_y*dt

            @. φ += ruidoDr + 5*(quorum/Nc)*dt

            @. x += vx
            @. y += vy

            x, y = periodic_bc(x,y,L)   
            
            
                if i % intervalo == 0
                    # Guardar datos cada 100 pasos de tiempo
                    x_data[Int64(i/intervalo),:]  .= x 
                    y_data[Int64(i/intervalo),:]  .= y
                    φ_data[Int64(i/intervalo),:]  .= φ
                    vx_data[Int64(i/intervalo),:] .= vx
                    vy_data[Int64(i/intervalo),:] .= vy
                end



               
                
                if i % 200 == 0

                    #actualizar la lista de vecinos cada ciertos pasos
                    vf, vc, dx, dy, r2 = vecinos_pbc(x,y,L,α)

                end
                
        
 


                #φ_old = φ
                #x_old = x
                #y_old = y 

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
        writedlm(joinpath(carpeta, "vx=$(round(v, digits=2)).csv"), vx_data, ',')
        writedlm(joinpath(carpeta, "vy=$(round(v, digits=2)).csv"), vy_data, ',')

    return
end



function correccion_lj(veci_fuerza, dx, dy, r2, radio, n_particulas)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    r_cut  = 2 * radio
    r2_cut = r_cut^2

    for i in 1:n_particulas
        for j in veci_fuerza[i]
            rij2 = r2[i, j]
            if rij2 == 0.0 || rij2 > r2_cut
                continue
            end

            rij = sqrt(rij2)

            f = soft_fuerza(rij, 75.0, r_cut)

            invr = 1 / rij
            fuerza_x[i] -= f * dx[i, j] * invr
            fuerza_y[i] -= f * dy[i, j] * invr
        end
    end

    return fuerza_x, fuerza_y
end



function soft_fuerza(r, epsilon, sigma)
    
    return epsilon * (1 - r/sigma)^(3/2)
end


function quorum_sensing(φ, veci_cono, dx, dy, r2,angulo1, α, Ro=3.0)
    N = length(φ)

    quorum = zeros(N)
    Nc = zeros(N)   

    cosφ = cos.(φ)
    sinφ = sin.(φ)
    cos_cono = cos(angulo1)

    r2_max = (α * Ro)^2

    for i in 1:N
        cφ = cosφ[i]
        sφ = sinφ[i]

        for j in veci_cono[i]
            rij2 = r2[i, j]
            if rij2 > r2_max
                continue
            end

            rij = sqrt(rij2)

            dxij = dx[i, j]
            dyij = dy[i, j]

            invr = 1 / rij
            rx = dxij * invr
            ry = dyij * invr

            if rx*cφ + ry*sφ >= cos_cono
                w = exp(-rij / Ro)
                ang = atan(ry, rx)

                quorum[i] += w * sin(ang - φ[i])
                Nc[i]     += w
            end
        end

        # Si no hay interaccion se fija a 1 para evitar division por cero
        if Nc[i] == 0.0
            Nc[i] = 1.0
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
    param_str = "N=$(n_particulas)-Dr=$(Dr)-θ=$(angle_str)-Ro=$(3*α)"
    
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

function generar_log(folder_path, v, n_pasos, n_particulas, L, angulo1, η, α, Dr, dt)
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
        println(file, "Tamaño caja: $L")
        println(file, "Empaquetamiento: $η")
        println(file, "Tamaño cono: $(3*α)")
    end
end



function packing_fraction(N, L, r=1)
    η = (π * N * r^2) / L^2
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


function vecinos_pbc(x, y,L,rq=4.0)
    n = length(x)
    invL = 1.0 / L

    veci_fuerza = Vector{Vector{Int}}(undef, n)
    veci_cono   = Vector{Vector{Int}}(undef, n)

    # --- relative vectors: r_ij = r_j - r_i ---
    dx = x' .- x
    dy = y' .- y

    # --- minimum image convention ---
    dx .-= L .* round.(dx .* invL)
    dy .-= L .* round.(dy .* invL)

    r2 = dx.^2 .+ dy.^2

    # --- cutoffs (squared) ---
    r_fuerza2 = (2 * 2^(1/6) + 1.0)^2
    r_quorum2 = (rq * 3 + 1.0)^2

    for i in 1:n
        veci_fuerza[i] = findall(j -> 0 < r2[i,j] ≤ r_fuerza2, 1:n)
        veci_cono[i]   = findall(j -> 0 < r2[i,j] ≤ r_quorum2, 1:n)
    end

    return veci_fuerza, veci_cono, dx, dy, r2
end


function ini_pbc(n_particulas::Int,L::Float64,α::Float64;radio::Float64 = 1.0,dt::Float64 = 0.01,pasos::Int = 3000,μ::Float64 = 1.0)

    # Posiciones aleatorias con overlap (en principio) y orientaciones
    φ = 2π .* rand(n_particulas)
    x = (rand(n_particulas) .- 0.5) .* L
    y = (rand(n_particulas) .- 0.5) .* L

    # Pasos de relajación
    for _ in 1:pasos

        # Lista de vecinos
        veci_fuerza, _,dx, dy, r2 = vecinos_pbc(x, y, L,α)
        # Fuerza
        Fx, Fy = correccion_lj(veci_fuerza, dx, dy, r2, radio, n_particulas)

        # Evolucionar el sistema
        @. x += μ * Fx * dt
        @. y += μ * Fy * dt

       
        
        #Aplicar PBC
        x, y = periodic_bc(x,y,L)   
    end

    return x, y, φ
end

