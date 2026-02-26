using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances, Base.Threads



dt = 1e-4  #Paso temporal

intervalo = 2000


function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, radio::Float64, angulo1::Float64,α::Float64,Dr::Float64,radio_p = 1,radio_b = 0.1)
   
     # Carpeta donde se guardara los datos de la simulacion
     carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α)
        η       = packing_fraction(n_particulas,radio)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, radio, angulo1,η,α,Dr,dt)
     # Guardar seed
        #Random.seed!(seed)
        n_updates = 0
        sqrtT = sqrt(2*Dr*dt)
        δ_fuerza = 0.5 * radio_p
        δ_cono   = 0.5 * radio_p
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        y_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vx_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vy_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        φ_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)

        vx, vy = zeros(n_particulas), zeros(n_particulas)
        
        x,y,φ = ini_circular(n_particulas,radio,α)
        x_ref, y_ref = copy(x), copy(y)
        vf, vc = vecinos(x,y,α)

        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y    = correccion_lj(vf,x,y,n_particulas)
            
            quorum, Nc  = quorum_sensing(φ,vc,x,y,angulo1,α)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)

            @. φ += ruidoDr + 5*(quorum/Nc)*dt

            @. vx = v*cos(φ)*dt + f_x*dt
            @. vy = v*sin(φ)*dt + f_y*dt

            @. x += vx
            @. y += vy

            barrera_circular!(x,y,radio,radio_p)
            
            if i % intervalo == 0
                # Guardar datos cada 100 pasos de tiempo
                x_data[Int64(i/intervalo),:]  .= x 
                y_data[Int64(i/intervalo),:]  .= y
                φ_data[Int64(i/intervalo),:]  .= φ
                vx_data[Int64(i/intervalo),:] .= vx
                vy_data[Int64(i/intervalo),:] .= vy
            end

            dmax2 = chequeo_lista(x, y, x_ref, y_ref)
                
            if dmax2 > δ_fuerza/2*δ_fuerza/2 || dmax2 > δ_cono/2*δ_cono/2
                vf, vc = vecinos(x, y, α)
            
                x_ref .= x
                y_ref .= y
                n_updates += 1
            end
           
        

        

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
        writedlm(joinpath(carpeta, "vx=$(round(v, digits=2)).csv"), vx_data, ',')
        writedlm(joinpath(carpeta, "vy=$(round(v, digits=2)).csv"), vy_data, ',')
        
    return n_updates
end




function correccion_lj(veci_fuerza, x, y,n_particulas,radio=1)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    r_cut  = 2^(1/6)*2 * radio
    r2_cut = r_cut^2

    
    Threads.@threads for i in 1:n_particulas

        xi, yi = x[i], y[i]
    
        fx = 0.0
        fy = 0.0
    
        for j in veci_fuerza[i]
            dx = x[j] - xi
            dy = y[j] - yi
            r2 = dx*dx + dy*dy
    
            if r2 == 0.0 || r2 > r2_cut
                continue
            end
    
            rij = sqrt(r2)
            f = lj_fuerza(rij, 0.5, 2*radio)
    
            fx += f * dx
            fy += f * dy
        end
    
        fuerza_x[i] = fx
        fuerza_y[i] = fy
    end

    return fuerza_x, fuerza_y
end

function correccion_soft(veci_fuerza, x, y, n_particulas, radio= 1.0)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    r_cut  = 2 * radio
    r2_cut = r_cut^2

    Threads.@threads for i in 1:n_particulas

        xi, yi = x[i], y[i]

        for j in veci_fuerza[i]
            dx = x[j] - xi
            dy = y[j] - yi
            r2 = dx^2 + dy^2

            if r2 == 0.0 || r2 > r2_cut
                continue
            end
    
            rij = sqrt(r2)

            f = soft_fuerza(rij, 75.0, r_cut)

            invr = 1 / rij
            fuerza_x[i] -= f * dx * invr
            fuerza_y[i] -= f * dy * invr
        end
    end    

    return fuerza_x, fuerza_y
end



function soft_fuerza(r, epsilon, sigma)
    
    return epsilon * (1 - r/sigma)^(3/2)
end

function lj_fuerza(distancia, epsilon, sigma)
    # Calculate the potential energy
    force = 24 * epsilon * ((sigma^6) / (distancia^8) - 2 * (sigma^12) / (distancia^14))
    return force
end




function quorum_sensing(φ, veci_cono, x, y, angulo1, α, Ro=3.0)
    N = length(φ)

    quorum = zeros(N)
    Nc = zeros(N)   

    cosφ = cos.(φ)
    sinφ = sin.(φ)
    cos_cono = cos(angulo1)

    r2_cut = (α * Ro)^2

    Threads.@threads for i in 1:N

        cφ = cosφ[i]
        sφ = sinφ[i]
        xi, yi = x[i], y[i]

        qi = 0.0
        Ni = 0.0

        for j in veci_cono[i]

            dx = x[j] - xi
            dy = y[j] - yi
            r2 = dx*dx + dy*dy

            if r2 == 0.0 || r2 > r2_cut
                continue
            end
    
            rij = sqrt(r2)
            

            invr = 1 / rij
            rx = dx * invr
            ry = dy * invr

            if rx*cφ + ry*sφ >= cos_cono
                w = exp(-rij / Ro)
                ang = atan(ry, rx)

                qi += w * sin(ang - φ[i])
                Ni     += w
            end
        end

        # Si no hay interaccion se fija a 1 para evitar division por cero
        if Ni == 0.0
            Ni = 1.0
        end

        quorum[i] = qi
        Nc[i]     = Ni


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

    # ---- Angulo bonito ----
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
        string(angulo1)
    end

    Ro = 3α

    # ---- Jerarquía de carpetas ----
    base_N   = joinpath(base_dir, "N=$(n_particulas)")
    base_ang = joinpath(base_N, "θ=$(angle_str)")
    base_Dr  = joinpath(base_ang, "Dr=$(Dr)")
    base_Ro  = joinpath(base_Dr, "Ro=$(Ro)")

    # Crea toda la jerarquía si no existe
    mkpath(base_Ro)

    # ---- Carpeta de simulación numerada ----
    sim_id = 1
    sim_dir = joinpath(base_Ro, "sim_$(lpad(sim_id, 3, '0'))")

    while isdir(sim_dir)
        sim_id += 1
        sim_dir = joinpath(base_Ro, "sim_$(lpad(sim_id, 3, '0'))")
    end

    mkdir(sim_dir)
    return sim_dir
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




function vecinos(x::Vector{Float64}, y::Vector{Float64}, rq=4.0, radio=1.0)
    skin = 1.0
    n = length(x)

    veci_fuerza = Vector{Vector{Int}}(undef, n)
    veci_cono   = Vector{Vector{Int}}(undef, n)

    r_fuerza = (2 * radio + skin)^2
    r_quorum = (rq * 1.5 * radio + skin)^2

    for i in 1:n
        xi, yi = x[i], y[i]

        vf = Int[]
        vc = Int[]
        sizehint!(vf, 32)
        sizehint!(vc, 64)

        for j in 1:n
            j == i && continue

            dx = x[j] - xi
            dy = y[j] - yi
            r2 = dx*dx + dy*dy

            if r2 ≤ r_fuerza
                push!(vf, j)
            end
            if r2 ≤ r_quorum
                push!(vc, j)
            end
        end

        veci_fuerza[i] = vf
        veci_cono[i]   = vc
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







function barrera_circular!(x, y, Rbox, a)
    Rmax = Rbox - a
    Rmax2 = Rmax^2
    N = length(x)

    Threads.@threads for i in 1:N
        r2 = x[i]^2 + y[i]^2
        if r2 > Rmax2
            scale = Rmax / sqrt(r2)
            x[i] *= scale
            y[i] *= scale
        end
    end
    return nothing
end



function ini_circular(n_particulas::Int, Rbox::Float64, α::Float64; radio::Float64 = 1.0, dt::Float64 = 0.01, pasos::Int = 6000, μ::Float64 = 1.0)
    x, y = zeros(n_particulas), zeros(n_particulas)
    # --- Posiciones aleatorias dentro del circulo
    for i in 1:n_particulas 
        θ = 2π * rand() 
        r = sqrt(rand()) * (Rbox - radio) 
        x[i] = r * cos(θ) 
        y[i] = r * sin(θ) 
    end

    # Orientacion aleatoria
    φ = 2π * rand(n_particulas)
    # Reflejar en caso de que sea necesario
    barrera_circular!(x, y, Rbox, radio)


    # --- Pasos de relajacion
    @showprogress "Acomodando las partículas..." for _ in 1:pasos

        veci_fuerza, _ = vecinos(x,y,α)
        Fx, Fy = correccion_soft(veci_fuerza, x, y, n_particulas)

        @. x += μ * Fx * dt
        @. y += μ * Fy * dt

        barrera_circular!(x, y, Rbox, radio)
    end

    return x, y, φ
end


function chequeo_lista(x, y, x_ref, y_ref)

    nt = Threads.nthreads()
    max_local = zeros(nt)

    Threads.@threads for i in eachindex(x)

        tid = Threads.threadid()

        dx = x[i] - x_ref[i]
        dy = y[i] - y_ref[i]
        d2 = dx*dx + dy*dy

        if d2 > max_local[tid]
            max_local[tid] = d2
        end
    end

    return maximum(max_local)
end




function distancia(x,y)

    dx = x' .- x      # dx[i,j] = x[j] - x[i]
    dy = y' .- y      # dy[i,j] = y[j] - y[i]
    r2 = dx.^2 .+ dy.^2


  return dx, dy, r2  
end