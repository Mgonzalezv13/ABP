using Random, Plots, DelimitedFiles, ProgressMeter

folder_path = "/home/mayron/ABP"


#Variables asociadas a la barrera circular
L = 40  #diametro 
centro_x = 0  # centro en  x
centro_y = 0  # centro en  y
radio = L/2   # radio
inicio_gap = 0  # inicio del gap
fin_gap = 0     # fin del gap 


Dt = 0.22   #Difusion Traslacional
Dr = 0.16   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10e-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)




function pos(v, Ω, n_pasos, n_particulas)
    x = zeros(n_pasos, n_particulas)
        y = similar(x)
        φ = similar(x)
        x[1, :] .= 0
        y[1, :] .= 0
        φ[1, :] .= 0
    @showprogress "Calculando las trayectorias..." for j in 1:n_particulas

        for i in 1:n_pasos-1
            ruidoDtx = sqrtD * randn()
            ruidoDty = sqrtD * randn()
            ruidoDr = sqrtT * randn()
            φ[i+1, j] = φ[i, j] + Ω * dt + ruidoDr
            x[i+1, j] = x[i, j] + v * cos(φ[i+1, j]) * dt + ruidoDtx
            y[i+1, j] = y[i, j] + v * sin(φ[i+1, j]) * dt + ruidoDty
        end

    end

    writedlm("/home/mayron/Datos/barrera_$v/pos_x.csv", x, ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y.csv", y, ',')

end