from turtle import color
from matplotlib import pyplot as plt
from matplotlib.collections import PathCollection
import matplotlib.animation as animation
from typing import NamedTuple, List
import json
import numpy as np


COLOUR_MAP = {
    0: "#31e600", # green
    1: "#a07047", # brown
    2: "#ff078e", # magenta
    3: "#0064ff", # blue
    4: "#f4e503" # yellow
}

COLOUR_MAP = {
    0: "#31e600", # green
    1: "#a07047", # brown
    2: "#ff078e", # magenta
    3: "#0064ff", # blue
    4: "#f4e503" # yellow
}

class Settings(NamedTuple):
    num_particles: int = 5000
    width: int = 250
    height: int = 250
    radius: int = 5
    alpha: int = 180
    beta: int = 17
    velocity: float = 0.67
    time_steps: int = 1000

class State(NamedTuple):
    coords: np.array
    colours: np.array

# global var for animation update
points = None
states: List[State] = []

def load_settings_json(json_fpath: str):

    with open(json_fpath, "r") as json_file:
        json_dict = json.load(json_file)

    return Settings(
        num_particles=json_dict["num_particles"],
        width=json_dict["width"],
        height=json_dict["height"],
        radius=json_dict["radius"],
        alpha=json_dict["alpha"],
        beta=json_dict["beta"],
        velocity=json_dict["velocity"],
        time_steps=json_dict["time_steps"]
    )


def read_binary_out_file(out_file_fp: str, time_steps: int, num_p: int) -> List[State]:
    """
    Load particle data from binary file.
    """
    global states

    print("Num p: ", num_p)

    # open binary data file
    with open(out_file_fp, "rb") as bin_file:
        fdata = np.fromfile(bin_file, dtype=np.double, count=-1)

    # read all states and their corresponding particle positions / colours  
    # states = []
    s_off = 0
    print(fdata.shape)
    for _ in range(time_steps):
        xs = fdata[s_off: s_off + num_p]
        ys = fdata[s_off + num_p: s_off + 2 * num_p]
        coords = np.column_stack((xs, ys))
        colours = np.array([COLOUR_MAP[int(x) % 5] for x in fdata[s_off + 2 * num_p: s_off + 3 * num_p]])
        s_off += 3 * num_p
        states.append(State(coords, colours))


def update(i: int):
    global points
    points.set_offsets(states[i].coords)
    points.set_facecolors(states[i].colours)
    points.set_edgecolors(states[i].colours)
    return points,


def run_animation():
    global points
    
    # initialize plot
    fig, ax = plt.subplots()
    ax.set_facecolor('black')
    points = ax.scatter(states[0].coords[:, 0], states[0].coords[:, 1], color=states[0].colours, s=5)

    # run animation
    ani = animation.FuncAnimation(fig, update, range(len(states)), interval=30)
    plt.show()
    ani.save('animation.gif', writer='Pillow', fps=30)


if __name__ == "__main__":

    # load settings json
    settings = load_settings_json("settings.json")

    # open and load log file
    read_binary_out_file("out_display.bin", settings.time_steps, settings.num_particles)

    # for s in states:
    # s = states[0]
    # for i in range(settings.num_particles):
    #     print(f"{s.xs[i]}, {s.ys[i]}, {s.colours[i]}")

    # run animation in matplotlib
    run_animation()
