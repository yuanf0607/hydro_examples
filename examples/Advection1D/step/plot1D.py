import glob
import matplotlib.pyplot as plt

def read_data(file_name):
    with open(file_name) as f:
        data = [tuple(map(float, line.split())) for line in f]
    return zip(*data)  # Transpose rows to columns


label = "Discontinuous step function"
def plot_data(file_names):
    max_y = float("-inf")
    plt.figure(figsize=(10, 5))
    for file_name in file_names:
        x, y = read_data(file_name)
        max_y = max(max_y, max(y))
        plt.plot(x, y, color='blue')
    plt.title(label)
    plt.xlabel("X-values")
    plt.ylabel("Y-values")
    plt.ylim(0, max_y * 1.2)
    #plt.legend()

    # Save the plot as a PNG file
    plot_name = label.lower().replace(" ", "_") + ".png"
    plt.savefig(plot_name)

file_names = glob.glob("*.asc")
plot_data(file_names)
