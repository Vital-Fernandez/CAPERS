import matplotlib.pyplot as plt
import numpy as np

# Example datasets
datasets = [np.random.rand(100) * i for i in range(1, 6)]  # 5 different datasets
responses = []  # Store (index, key pressed)

fig, ax = plt.subplots()
index = [0]  # mutable container to store current dataset index

line, = ax.plot(datasets[0])
ax.set_title(f"Dataset 0")

def on_key(event):
    key = event.key
    i = index[0]

    # Record response
    responses.append((i, key))

    print(f"Pressed: {key} for dataset {i}")

    # Advance to next dataset
    i += 1
    if i >= len(datasets):
        plt.close(fig)
        print("Finished all datasets.")
        print("Responses:", responses)
        return

    # Update plot
    index[0] = i
    line.set_ydata(datasets[i])
    ax.set_title(f"Dataset {i}")
    ax.relim()
    ax.autoscale_view()
    fig.canvas.draw()

fig.canvas.mpl_connect('key_press_event', on_key)
plt.show()