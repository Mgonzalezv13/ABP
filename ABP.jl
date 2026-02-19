using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-4  #Paso temporal


intervalo = 2000

function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, L::Float64, angulo1::Float64,α,Dr,radio_p = 1)
   
     # Carpeta donde se guardara los datos de la simulacion
        carpeta = carpeta_simulacion("/home/mayron/Datos", angulo1,n_particulas,Dr,α,L)
        η       = packing_fraction(n_particulas,L)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, L, angulo1,η,α,Dr,dt)

        n_updates = 0
        sqrtT = sqrt(2*Dr*dt) #esta cantidad se mantiene fija
        δ_fuerza = 0.5 * radio_p
        δ_cono   = 0.5 * radio_p

        x_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        y_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vx_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        vy_data  = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)
        φ_data   = Matrix{Float64}(undef, Int64(n_pasos/intervalo), n_particulas)

        vx, vy = zeros(n_particulas), zeros(n_particulas)
        x, y, φ = ini_con_pbc(n_particulas,L,α)
        x_ref, y_ref = copy(x), copy(y)
        vf,vc = vecinos(x, y, L, α)

        @showprogress "Calculando..." for i in 2:n_pasos

            f_x, f_y = correccion_soft(vf,x, y,n_particulas,L)
            quorum, Nc = quorum_sensing(φ,vc,x,y,angulo1,α,L)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)

            @. φ += ruidoDr + 5*(quorum/Nc)*dt

            @. vx = v*cos(φ)*dt + f_x*dt
            @. vy = v*sin(φ)*dt + f_y*dt

            @. x += vx
            @. y += vy

            pbc!(x,y,L)
            
            
                if i % intervalo == 0
                    
                    x_data[Int64(i/intervalo),:]  .= x 
                    y_data[Int64(i/intervalo),:]  .= y
                    φ_data[Int64(i/intervalo),:]  .= φ
                    vx_data[Int64(i/intervalo),:] .= vx
                    vy_data[Int64(i/intervalo),:] .= vy
                end

                dmax2 = chequeo_lista(x, y, x_ref, y_ref,L)
                
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

function correccion_lj(veci_fuerza, x, y,n_particulas,L,radio=1)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    r_cut  = 2^(1/6)*2 * radio
    r2_cut = r_cut^2

    
    for i in 1:n_particulas

        xi, yi = x[i], y[i]
    
        fx = 0.0
        fy = 0.0
    
        for j in veci_fuerza[i]
            dx = x[j] - xi
            dy = y[j] - yi
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)
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

function correccion_soft(veci_fuerza, x, y, n_particulas, L, radio= 1.0)
    fuerza_x = zeros(n_particulas)
    fuerza_y = zeros(n_particulas)

    r_cut  = 2 * radio
    r2_cut = r_cut^2

    for i in 1:n_particulas

        xi, yi = x[i], y[i]

        for j in veci_fuerza[i]
            dx = x[j] - xi
            dy = y[j] - yi
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)

            r2 = dx*dx + dy*dy

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


function quorum_sensing(φ, veci_cono, x, y, angulo1, α,L, Ro=3.0)
    N = length(φ)

    quorum = zeros(N)
    Nc = zeros(N)   

    cosφ = cos.(φ)
    sinφ = sin.(φ)
    cos_cono = cos(angulo1)

    r2_cut = (α * Ro)^2

    for i in 1:N

        cφ = cosφ[i]
        sφ = sinφ[i]
        xi, yi = x[i], y[i]

        qi = 0.0
        Ni = 0.0

        for j in veci_cono[i]

            dx = x[j] - xi
            dy = y[j] - yi
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)

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


function vecinos(x::Vector{Float64}, y::Vector{Float64}, L , rq=4.0, radio=1.0)
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
            dx -= L * round(dx / L)
            dy -= L * round(dy / L)
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

function ini_con_pbc(n_particulas::Int,L::Float64,α::Float64;radio::Float64 = 1.0,dt::Float64 = 0.001,pasos::Int = 6000,μ::Float64 = 1.0)

    # --- Posiciones aleatorias en la caja
    x = L * (rand(n_particulas) .- 0.5)
    y = L * (rand(n_particulas) .- 0.5)

    # --- Orientaciones aleatorias
    φ = 2π * rand(n_particulas)

    # Asegurar PBC inicial
    pbc!(x, y, L)

    # --- Relajación
    @showprogress "Acomodando las partículas..." for _ in 1:pasos

        veci_fuerza, _ = vecinos(x, y, L,α)   # ← debe usar MIC
        Fx, Fy = correccion_soft(veci_fuerza, x, y, n_particulas,L)

        @. x += μ * Fx * dt
        @. y += μ * Fy * dt

        pbc!(x, y, L)
    end

    return x, y, φ
end


function pbc!(x, y, L)
    @. x = mod(x + L/2, L) - L/2
    @. y = mod(y + L/2, L) - L/2
    return nothing
end

function chequeo_lista(x, y, x_ref, y_ref,L)
    maxd2 = 0.0
    @inbounds for i in eachindex(x)
        dx = x[i] - x_ref[i]
        dy = y[i] - y_ref[i]
        dx -= L * round(dx / L)
        dy -= L * round(dy / L)
        d2 = dx*dx + dy*dy
        if d2 > maxd2
            maxd2 = d2
        end
    end
    return maxd2
end



function carpeta_simulacion(base_dir, angulo1, n_particulas, Dr, α,L)

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
    base_N   = joinpath(base_dir, "PBC_N=$(n_particulas)")
    base_L   = joinpath(base_N, "Tamaño_Caja=$(L)")
    base_ang = joinpath(base_L, "θ=$(angle_str)")
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