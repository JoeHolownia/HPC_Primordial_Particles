from turtle import color
from matplotlib import pyplot as plt
from matplotlib.collections import PathCollection
import matplotlib.animation as animation
from typing import NamedTuple, List
import json
import numpy as np
import os


COLOUR_MAP = {
    0: "#31e600", # green
    1: "#a07047", # brown
    2: "#ff078e", # magenta
    3: "#0064ff", # blue
    4: "#f4e503" # yellow
}

class Settings(NamedTuple):
    num_particles: int
    width: int = 250
    height: int = 250
    density: int = 8
    alpha: int = 180
    beta: int = 17
    gamma: float = 0.134
    time_steps: int = 1000

class State:

    def __init__(self, time_step):
        self.time_step = time_step
        self.coords = []
        self.colours = []

    def append_data(self, new_coords: list, new_colours: np.array):
        self.coords.extend(new_coords)
        self.colours.extend(new_colours)

    def get_coords_np(self):
        return np.array(self.coords)

    def get_colours_np(self):
        return np.array(self.colours)

    def get_time_step(self):
        return self.time_step


# global vars for animation update
points = None
states: List[State] = []
text: str = "Step 0"
fig, ax = plt.subplots()

def load_settings_json(json_fpath: str):
    """
    Loads the simulation settings from json into a settings struct.
    """

    with open(json_fpath, "r") as json_file:
        json_dict = json.load(json_file)

    return Settings(
        num_particles=json_dict["num_particles"],
        width=json_dict["width"],
        height=json_dict["height"],
        density=json_dict["density"],
        alpha=json_dict["alpha"],
        beta=json_dict["beta"],
        gamma=json_dict["gamma"],
        time_steps=json_dict["time_steps"]
    )


def read_binary_out_file(bin_file_dir: str, time_steps: int, num_p: int) -> List[State]:
    """
    Load particle data from binary file.
    """
    global states
    states = [State(i) for i in range(time_steps)]

    # loop through binary files and read
    for out_fname in sorted(os.listdir(bin_file_dir)):

        # skip non-binary files
        if out_fname.split(".")[1] != "bin":
            continue

        # DEBUG!!!
        # if out_fname != "out_display_3.bin":
        #     continue

        # open binary data file
        out_file_fp = bin_file_dir + "/" + out_fname
        with open(out_file_fp, "rb") as bin_file:
            fdata = np.fromfile(bin_file, dtype=np.single, count=-1)

        # read all states and their corresponding particle positions / colours  
        s_off = 0
        for i in range(time_steps):
            num_p = int(fdata[s_off])
            part_off = s_off + 1
            xs = fdata[part_off: part_off + num_p]
            ys = fdata[part_off + num_p: part_off + 2 * num_p]
            coords = [(xs[i], ys[i]) for i in range(len(xs))]
            colours = [COLOUR_MAP[int(x) % 5] for x in fdata[part_off + 2 * num_p: part_off + 3 * num_p]]
            s_off += 1 + 3 * num_p
            states[i].append_data(coords, colours)

def update(i: int):
    """
    Update function for matplotlib animation.
    """
    global points
    global text
    points.set_offsets(states[i].get_coords_np())
    points.set_facecolors(states[i].get_colours_np())
    points.set_edgecolors(states[i].get_colours_np())
    text.set_text(f"Step: {states[i].get_time_step()}")
    return points,


def run_animation(settings: Settings):
    """
    Function which sets up plot and runs animation with matplotlib.
    """

    global points
    global text
    
    # initialize plot
    ax.set_facecolor('black')
    ax.set_xlim((0, settings.width))
    ax.set_ylim((0, settings.height))
    points = ax.scatter(states[0].get_coords_np()[:, 0], states[0].get_coords_np()[:, 1], color=states[0].get_colours_np(), s=5)
    text = ax.text(0.22, 0.96, "Step: 0", color='white', fontsize=14, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes)

    # save png of initial state
    fig.savefig("../results/init_state.png", format='png')

    # run animation
    ani = animation.FuncAnimation(fig, update, range(len(states)), interval=60)
    plt.show()

    print("Saving animation to gif...")
    ani.save('../results/animation.gif', writer='Pillow', fps=30)
    print("Saved to animation to results/animation.gif")


if __name__ == "__main__":

    # load settings json
    settings = load_settings_json("settings.json")

    # open and load display files
    read_binary_out_file("../results", settings.time_steps + 1, settings.num_particles)

    # run animation in matplotlib
    run_animation(settings)
