import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker

def read_files(data_path="./data"):
    # Read files from data path, expected on the form 'e_xy_X.csv' 
    # with 'e_' or 'r_', then 'xy_' or 'xz_', then 'X' (number).
    e_xy = []
    e_xz = []
    r_xy = []
    r_xz = []
    for filename in os.listdir(data_path):
        arr = pd.read_csv(f"{data_path}/{filename}", header=None).values
        if re.match(r"^(e_xy_)(\d)+(.csv)$", filename):
            e_xy.append(arr)
        elif re.match(r"^(e_xz_)(\d)+(.csv)$", filename):
            e_xz.append(arr)
        elif re.match(r"^(r_xy_)(\d)+(.csv)$", filename):
            r_xy.append(arr)
        elif re.match(r"^(r_xz_)(\d)+(.csv)$", filename):
            r_xz.append(arr)
        else:
            raise NameError(f"Invalid filename '{filename}' found in data folder.")
    
    return [[e_xy, e_xz], [r_xy, r_xz]]


def animate(arrays):
    # Animate e and r as two subplots of imshow (image representation of a 2D mesh)
    fig, ax = plt.subplots(2, 2)
    plt.rc('axes', titlesize=8)
    e_xy = arrays[0][0]
    frames = len(e_xy)
    dim = e_xy[0].shape[0]
    titles = [["e front (xy)", "e right (xz)"], ["r front (xy)", "r right (xz)"]]
    for i in range(2):
        for j in range(2): 
            ax[i, j].set_title(titles[i][j])
            im = ax[i, j].imshow(arrays[i][j][0], vmin=0.0, vmax=1.0)
            ax[i, j].xaxis.set_major_locator(ticker.NullLocator())
            ax[i, j].yaxis.set_major_locator(ticker.NullLocator())
            fig.colorbar(im, ax=ax[i][j])

    # For matplotlibs FuncAnimation to call for each new frame
    def new_frame(num):
        for i in range(2):
            for j in range(2):
                ax[i, j].imshow(arrays[i][j][num], vmin=0.0, vmax=1.0)

    ani = animation.FuncAnimation(
        fig, new_frame, frames, interval=200, blit=False
    )

    ani.save(f"./result.gif",  fps=5,  writer="imagemagick")


if __name__ == "__main__":
    arrays = read_files()
    animate(arrays)
