using Random, Serialization, Plots, DelimitedFiles

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 1000000 #se bajo el numero de pasos temporales, debido a que con muchos colapsaba el grafico
n_particulas = 5

#Aca se definen las constantes 
# Define box parameters
L = 100  # length of box side
x_left = -L/2  # barrera izquierda
x_right = L/2  # barrera derecha
y_bottom = -L/2  # barrera inferior
y_top = L/2  # barrera superior

v  = 20
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
            #msd[i] =(x[i+1] .- x[1]).^2 + (y[i+1] .- y[1]).^2
            #msd[1000000] = (x[1000000] .- x[1]).^2 + (y[1000000] .- y[1]).^2
            # Add reflective boundary conditions for x-coordinate
        #if x[i+1] < x_left
           # x[i+1] = x_left + abs(x[i+1] - x_left)  # Reflect particle off left wall
           # φ[i+1] = pi - φ[i+1]  # Reflect angle
       # elseif x[i+1] > x_right
          #  x[i+1] = x_right - abs(x[i+1] - x_right)  # Reflect particle off right wall
          #  φ[i+1] = pi - φ[i+1]  # Reflect angle
       # end
        
        # Add reflective boundary conditions for y-coordinate
        #if y[i+1] < y_bottom
           # y[i+1] = y_bottom + abs(y[i+1] - y_bottom)  # Reflect particle off bottom wall
            #φ[i+1] = -φ[i+1]  # Reflect angle
       # elseif y[i+1] > y_top
          #  y[i+1] = y_top - abs(y[i+1] - y_top)  # Reflect particle off top wall
          #  φ[i+1] = -φ[i+1]  # Reflect angle
        end
        #end

    #Guarda la posicion en "traj$j_dat" y guarda el msd en "msd$j_dat" de cada particula con j el numero de la particula 
    file_path = joinpath(folder_path, "traj_$j.dat")
    writedlm(file_path, [x y])
    scatter(x,y, markersize=1, legend=false)
    savefig("traj_$j.png")
end









