import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_files(data_path="./data"):
    # Read files from data path
    # Filenames are expected on the form 'eX.csv' or 'rX.csv', 
    # where X is an arbitrary number
    e_list = []
    r_list = []
    for filename in os.listdir(data_path):
        arr = pd.read_csv(f"{data_path}/{filename}").values
        if re.match(r"^(e)(\d)+(.csv)$", filename):
            e_list.append(arr)
        elif re.match(r"^(r)(\d)+(.csv)$", filename):
            r_list.append(arr)
        else:
            raise NameError(f"Invalid filename: '{filename}' found in data folder.")
    
    return e_list, r_list


def animate(e_list, r_list):
    # Animate e and r as two subplots of imshow (image representation of a 2D mesh)
    fig, ax = plt.subplots(1, 2)

    ax[0].set_title("e mesh")
    ax[1].set_title("r mesh")
    im_e = ax[0].imshow(e_list[0])
    im_r = ax[1].imshow(r_list[0])
    fig.colorbar(im_e, ax=ax[0])
    fig.colorbar(im_r, ax=ax[1])

    # For matplotlibs FuncAnimation to call for each new frame
    def new_frame(num):
        ax[0].imshow(e_list[num])
        ax[1].imshow(r_list[num])

    ani = animation.FuncAnimation(
        fig, new_frame, len(e_list), interval=200, blit=False
    )

    ani.save("./result.gif",  fps=5,  writer="imagemagick")


if __name__ == "__main__":
    e_list, r_list = read_files()
    animate(e_list, r_list)

"""
def animate(e_list, r_list):
    # Animate e and r as two subplots of imshow (image representation of a 2D mesh)
    fig = plt.figure()
    ax_e = fig.add_subplot(121)
    ax_e.set_title("e mesh")
    im_e = ax_e.imshow(e_list[0])
    ax_r = fig.add_subplot(122)
    ax_r.set_title("r mesh")
    im_r = ax_r.imshow(r_list[0])
    fig.colorbar(im_r, ax=ax_r)

    # For matplotlibs FuncAnimation to call for each new frame
    def new_frame(num):
        im_e = ax_e.imshow(e_list[num])
        im_r = ax_r.imshow(r_list[num])

    ani = animation.FuncAnimation(
        fig, new_frame, len(e_list), interval=200, blit=False
    )

    ani.save("./result.gif",  fps=5,  writer="imagemagick")
"""