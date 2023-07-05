using Random, Plots, DelimitedFiles, ProgressMeter, Distributions

folder_path = "/home/mayron/ABP"

#aca se define el numero de componentes en la direccion x e y

L = 40  #diametro del circulo
centro_x = 0  # centro del circulo en x
centro_y = 0  # centro del circulo en y
radio = L/2  # radio del circulo


Dt = 0.22   #Difusion Traslacional
Dr = 0.16   #Difusion Rotacional
Ω  = 0.0    #Constante de quiralidad   
dt = 10^-3  #Paso temporal
sqrtD = sqrt(2*Dt*dt) #esto corresponde a √(2*Dt*dt)
sqrtT = sqrt(2*Dr*dt) #esto corresponde a √(2*Dr*dt)
uniform_dist = Uniform(0, 2π)

function barrera(v,Ω, n_pasos, n_particulas)
    pos_x = zeros(n_pasos,n_particulas)
    pos_y = similar(pos_x)
    @showprogress "Calculando trayectorias " for j in 1:n_particulas
        #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
            x   = zeros(n_pasos)
            y   = similar(x)
            φ   = similar(x)
            x[1] = 0
            y[1] = 0
            φ[1] = rand(uniform_dist)
            for i in 1:n_pasos-1
                
                ruidoDtx  = sqrtD*randn()
                
                ruidoDty  = sqrtD*randn() #Se definen ruidos diferentes para X e Y
                
                ruidoDr  = sqrtT*randn()
                
                φ[i+1] = φ[i] + Ω*dt +  ruidoDr
                
                x[i+1] = x[i] + v*cos(φ[i])*dt + ruidoDtx
                
                y[i+1] = y[i] + v*sin(φ[i])*dt +  ruidoDty
    
                # Angulo de la particula con respecto al centro del circulo
                dx = x[i+1] - centro_x
                dy = y[i+1] - centro_y
                angulo = atan(dy, dx)

                # Se ajusta el angulo entre -pi y pi
                angulo = angulo % (2*pi)
                if angulo > pi
                    angulo -= 2*pi
                end

                #Verificar si la particula esta fuera de la barrera
                distancia = sqrt((x[i+1] - centro_x)^2 + (y[i+1] - centro_y)^2)
                if distancia >= radio
                    # Reflejar la particula en la barrera
                    normal_x = (x[i+1] - centro_x) / distancia
                    normal_y = (y[i+1] - centro_y) / distancia
                    
                    
                    producto_punto = (v * cos(φ[i]) * normal_x) + (v * sin(φ[i]) * normal_y)
                    
                    # Calculate the reflected direction
                    reflect_x = 2 * producto_punto * normal_x - cos(φ[i])
                    reflect_y = 2 * producto_punto * normal_y - sin(φ[i])
                    
                    # Normalize the reflected direction
                    reflect_norm = sqrt(reflect_x^2 + reflect_y^2)
                    reflect_x /= reflect_norm
                    reflect_y /= reflect_norm
                    
                    # Update the position after reflection
                    x[i+1] = centro_x + reflect_x * radio
                    y[i+1] = centro_y + reflect_y * radio
                    
                    # Calculate the new angle
                    φ[i+1] = atan(reflect_y, reflect_x)
                end
            
            pos_x[:,j] = x
            pos_y[:,j] = y
            end
                

    end
    writedlm("/home/mayron/Datos/barrera_$v/pos_x_v=00$v.csv",pos_x , ',')
    writedlm("/home/mayron/Datos/barrera_$v/pos_y_v=00$v.csv",pos_y , ',')
end


function msd(v,Ω, n_pasos, n_particulas)
    msd_total = zeros(n_pasos,n_particulas)
    @showprogress "Calculando las trayectorias..." for j in 1:n_particulas
            #Aca se definen vectores "vacios" para almacenar las posiciones en x e y de cada particula 
                x   = zeros(n_pasos)
                y   = similar(x)
                φ   = similar(x)
                msd = similar(x)
                x[1]    = randn()
                y[1]    = randn()
                φ[1]    = rand(Uniform(0, 2π))
                        for i in 1:n_pasos-1
                            
                            ruidoDtx  = sqrtD*randn()
                            
                            ruidoDty  = sqrtD*randn() #Se definen ruidos diferentes para X e Y
                            
                            ruidoDr  = sqrtT*randn()
                            
                            φ[i+1] = φ[i] + Ω*dt +  ruidoDr
                            
                            x[i+1] = x[i] + v*cos(φ[i])*dt + ruidoDtx
                            
                            y[i+1] = y[i] + v*sin(φ[i])*dt +  ruidoDty

                            msd[i+1] = ( x[i+1]- x[1] )^2 + ( y[i+1] - y[1] )^2
                        end
                        msd[1]  = ( x[2]- x[1] )^2 + ( y[2] - y[1] )^2
                        msd_total[:,j] = msd
        end
        msd_prom = mean(msd_total, dims=2)
        writedlm("/home/mayron/Datos/datos2.0_$v/msd_00$v.csv",msd_prom, ',')
end










