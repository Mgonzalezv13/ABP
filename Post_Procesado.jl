using DelimitedFiles, Plots

folder_path = "/home/mayron/ABP/"

n_trayectorias = 10000
n_particulas = 50
dt = 10e-3

# Define the circular barrier parameters
L = 200
centro_x = 0
centro_y = 0
radio = L/2
inicio_gap = 0
fin_gap = 0

# Define the fill range for the circular barrier
theta_fill = LinRange(0, 2π, 100)
circ_x_fill = centro_x .+ radio * cos.(theta_fill)
circ_y_fill = centro_y .+ radio * sin.(theta_fill)

# Define the line range for the gap in the circular barrier
theta_gap = LinRange(inicio_gap, fin_gap, 100)
gap_x = centro_x .+ radio * cos.(theta_gap)
gap_y = centro_y .+ radio * sin.(theta_gap)


for j in 1:n_particulas
        # Loop over the particles
        Plots.closeall()

        # Preallocate arrays for position data
        x = zeros(n_trayectorias)
        y = similar(x)

        # Load position data from file
        pos = readdlm("traj_$j.dat")
        x = pos[:, 1]
        y = pos[:, 2]

        # Plot the circular barrier
        plot(circ_x_fill, circ_y_fill, fillrange = 0, fillalpha = 0.3, color = :gray, aspect_ratio = :equal, legend = false)

        # Plot the gap in the circular barrier
       # plot!(gap_x, gap_y, color = :white, aspect_ratio = :equal, linewidth = 15, legend = false)

        # Plot the particle trajectory
        scatter!(x, y, markersize = 0.2, legend = false)

        # Save the plot as an SVG file
        savefig(joinpath(folder_path, "traj_$j.svg"))
end