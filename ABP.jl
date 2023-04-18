using Distributions, LinearAlgebra, GLM, DataFrames, Printf, CSV
using Plots, Plots.PlotMeasures, StatsPlots
import Random


#aca se define el numero de componentes en la direccion x e y
n_trayectorias = 10000
n_particulas = 10

#Aca se definen las constantes 
L = 40
v  = fill(10.0, n_particulas)
Dt = fill(0.08, n_particulas)
Dr = fill(0.16, n_particulas)
Ω  = fill(0.0, n_particulas)
dt = 10e-3

#Aca se definen vectores "vacios
x  = zeros(n_trayectorias, n_particulas)
y  = similar(x)
φ  = similar(x)
Δt = zeros(n_trayectorias)
msd = similar(Δt)
x[1,:] .= 0
y[1,:] .= 0
φ[1,:] .= 0
Δt[1]   = 0
for i in 1:n_trayectorias-1
    for j in 1:n_particulas
        x[i+1,j] = x[i,j] + v[j]*cos(φ[i,j])*dt + sqrt(2*Dt[j]*Dr[j]*dt)*randn()
        y[i+1,j] = y[i,j] + v[j]*sin(φ[i,j])*dt + sqrt(2*Dt[j]*Dr[j]*dt)*randn()
        φ[i+1,j] = φ[i,j] + Ω[j]*dt + sqrt(2*Dt[j]*Dr[j]*dt)*randn()
        Δt[i+1]  =+ Δt[i] .+ dt*i
    end
end



# create a matrix and header

header = ["Particula_1", "Particula_2", "Particula_3", "Particula_4", "Particula_5", "Particula_6", "Particula_7", "Particula_8", "Particula_9", "Particula_10"]

# convert the matrix to a DataFrame
df = DataFrame(x, header)

# define the path and filename for the CSV file
filename = "posicion_y.csv"

# write the DataFrame to the CSV file
CSV.write(filename, df)

header = ["Particula_1", "Particula_2", "Particula_3", "Particula_4", "Particula_5", "Particula_6", "Particula_7", "Particula_8", "Particula_9", "Particula_10"]

# convert the matrix to a DataFrame
df = DataFrame(y, header)

# define the path and filename for the CSV file
filename = "posicion_y.csv"

# write the DataFrame to the CSV file
CSV.write(filename, df)





plot(x,y)

#for j in 1:n_particulas
 #   msd[j] = mean((x[:,j] .- x[1,j]).^2 + (y[:,j] .- y[1,j]).^2)
#end



#Grafico del msd



