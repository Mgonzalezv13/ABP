using  DelimitedFiles, LinearAlgebra, Printf, Dates, Clustering, Statistics, ProgressMeter, Random, Distances



dt = 1e-3  #Paso temporal




function vc(v::Int64, n_pasos::Int64, n_particulas::Int64, radio::Int64, angulo1::Float64,α,Dr,radio_p = 1)
   
     # Carpeta donde se guardara los datos de la simulacion
        carpeta = carpeta_simulacion("/home/mayron/Datos")
        η       = packing_fraction(n_particulas,radio)

     # Archivo log con los parámetros de la simulacion
        generar_log(carpeta, v, n_pasos, n_particulas, radio, angulo1,η,α,Dr)
     # Guardar seed
        #Random.seed!(seed)
   
        sqrtT = sqrt(2*Dr*dt)
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x_data  = Vector{Float64}[]
        y_data  = Vector{Float64}[]
        φ_data  = Vector{Float64}[]
        φ_old = rand(0:2pi,n_particulas)
        x_old , y_old = condicion_inicial(n_particulas,radio)
        vf,vc = vecinos(x_old,y_old,α)
        rastro_x,rastro_y = rastro(x_old,y_old)
        @showprogress "Calculando..." for i in 2:n_pasos

          
            f_x, f_y = correccion_lj(x_old, y_old, vf,radio_p,n_particulas)
            quorum, Nc = quorum_sensing(x_old, y_old, n_particulas,φ_old, angulo1,vc,α)
            

            
            ruidoDr  = sqrtT * randn(n_particulas)
            
            φ = φ_old + 5*(quorum./Nc)*dt + ruidoDr 
            
            x = x_old + v*cos.(φ_old)*dt  + f_x*dt
            
            y = y_old + v*sin.(φ_old)*dt  + f_y*dt


            
            
                if i % 100 == 0
                    # Guardar datos cada 100 pasos de tiempo
                    push!(x_data, x)
                    push!(y_data, y)
                    push!(φ_data, φ)
                end



               
                
                if i % 150 == 0

                    #actualizar la lista de vecinos cada ciertos pasos
                    vf,vc = vecinos(x,y,4)

                end
                
        

                
            x, y = reflective_bc(x_old, y_old, x, y, radio)


                φ_old = φ
                x_old = x
                y_old = y 

            end  
        # Guardar las posiciones en los archivos
        writedlm(joinpath(carpeta, "pos_x_v=$(round(v, digits=2)).csv"), x_data, ',')
        writedlm(joinpath(carpeta, "pos_y_v=$(round(v, digits=2)).csv"), y_data, ',')
        writedlm(joinpath(carpeta, "phi_v=$(round(v, digits=2)).csv"), φ_data, ',')
    return
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
    return fuerza_x, fuerza_y
end

function lj_fuerza(distancia, epsilon, sigma)
        # Calculate the potential energy
        force = 24 * epsilon * ((sigma^6) / (distancia^8) - 2 * (sigma^12) / (distancia^14))
        return force
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

function generar_log(folder_path, v, n_pasos, n_particulas, radio, angulo1, η,α,Dr)
    # Formatear la velocidad para incluirlo como string
    v_str = @sprintf("%.2f", v)
    R0 = 3*α
    # Escribir los parámetros en el archivo log
    log_filename = joinpath(folder_path, "log.txt")
    open(log_filename, "w") do file
        println(file, "Esta simulación se realizó el: $(Dates.now())")
        println(file, "Parámetros de la simulación:")
        println(file, "v: $v")
        println(file, "Número de iteraciones: $n_pasos")
        println(file, "Número de partículas: $n_particulas")
        println(file, "Radio de la barrera: $radio")
        println(file, "Angulo de apertura del cono de visión: $angulo1")
        println(file, "Fracción de empaquetamiento: $η")
        println(file, "Tamaño del cono: $R0")
        println(file, "Ruido Rotacional: $Dr")


        println(file, "--------------------------------------")
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