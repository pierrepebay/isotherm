import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 50  # Size of the temperature grid
num_iterations = 100  # Number of iterations to visualize

def read_temperature_data(iteration):
    filename = f"cmake-build-debug/output/temperature_data_{iteration}.txt"
    data = np.loadtxt(filename, delimiter=",")
    return data

# Set up the figure and axis for animation
fig, ax = plt.subplots()
cax = ax.imshow(np.zeros((N, N)), vmin=0, vmax=1000, cmap='plasma')
fig.colorbar(cax)

def animate(i):
    data = read_temperature_data(i)
    cax.set_array(data)
    ax.set_title(f"Iteration {i}")

# Create the animation
ani = animation.FuncAnimation(fig, animate, frames=num_iterations, interval=50)

# To save the animation, uncomment the following line
# ani.save('temperature_evolution.mp4', writer='ffmpeg')

plt.show()
