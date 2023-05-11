using Random, Serialization, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 100000 #se bajo el numero de pasos temporales, debido a que con muchos colapsaba el grafico
n_particulas = 5

L = 200  #diametro del circulo
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
radio = L/2  # radio del circulo
inicio_gap = 0  # inicio del agujero en la  barrera circular (en radianes)
fin_gap = pi/10 # fin del agujero en la barrera circular (en radianes)

v  = 5.5
Dt = 0.08
Dr = 0.16
Ω  = 0.0
dt = 10e-3
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2Dt*Dr*dt)
sqrtT = sqrt(2*Dt*dt)


 for j in 1:n_particulas
    #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
        x   = zeros(n_trayectorias)
        y   = similar(x)
        φ   = similar(x)
        x[1]    = 0
        y[1]    = 0
        φ[1]    = 0
        for i in 1:n_trayectorias-1
            ruidoDt  = sqrtD*randn() 
            ruidoDr  = sqrtT*randn()
            x[i+1] = x[i] .+ v*cos(φ[i])*dt .+ ruidoDt
            y[i+1] = y[i] .+ v*sin(φ[i])*dt .+  ruidoDt
            φ[i+1] = φ[i] .+ Ω*dt .+  ruidoDr

            # Angulo de la particula con respecto al centro del circulo
            dx = x[i+1] - centro_x
            dy = y[i+1] - centro_y
            angulo = atan(dy, dx)

            # Adjust angle to be within the range of -pi to pi
            angulo = angulo % (2*pi)
            if angulo > pi
                angulo -= 2*pi
            end

            # Check if the particle is within the gap region
            if angulo >= inicio_gap && angulo <= fin_gap
                continue  # Skip reflection if the particle is within the gap region
            end
                # Mientras la particula este dentro del circulo y fuera del angulo del gap
            distancia = sqrt(dx^2 + dy^2)
                
            if distancia > radio
                    # Reflect the particle off the circular wall
                    normal_x = dx / distancia  # x-component of outward normal vector
                    normal_y = dy / distancia  # y-component of outward normal vector
                    
                    # Reflect the particle's position
                    reflejo_x = centro_x + normal_x * radio
                    reflejo_y = centro_y + normal_y * radio
                    delta_x = reflejo_x - x[i+1]
                    delta_y = reflejo_y - y[i+1]
                    x[i+1] = reflejo_x + delta_x
                    y[i+1] = reflejo_y + delta_y
                    
                    # Reflect the angle
                    angle = atan(delta_y, delta_x)
                    φ[i+1] = angulo + pi + (angulo - φ[i])
            end
        end

    #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
    file_path = joinpath(folder_path, "traj_$j.dat")
    writedlm(file_path, [x y])
    # rellenar la zona exterior a la barrera, por contraste
    theta = LinRange(0, 2π, 100)
    circ_x = centro_x .+ radio * cos.(theta)
    circ_y = centro_y .+ radio * sin.(theta)
    plot!(circ_x, circ_y, fillrange = 0, fillalpha = 0.3, color = :gray, aspect_ratio = :equal, legend =false)
    theta = LinRange(inicio_gap, fin_gap, 100)
    gap_x = centro_x .+ radio * cos.(theta)
    gap_y = centro_y .+ radio * sin.(theta)
    plot!(gap_x, gap_y, color = :white, aspect_ratio = :equal, linewidth= 15, legend =false)
    #ploteo de la trayectoria de la particula
    scatter!(x,y, markersize=0.2, legend=false)
    savefig("traj_$j.pdf")
end









