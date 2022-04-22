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
    info = {"times": []}
    for filename in os.listdir(data_path):
        # Either read info file
        if filename == "info.txt":
            infile = open(f"{data_path}/{filename}")
            for line in infile.readlines():
                key, val = line.split("=")
                if re.match(r"^(t)(\d)+$", key):
                    info["times"].append(float(val))
                else:
                    info[key] = float(val)
            infile.close()
        else: # Or read csv files
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
    
    return [[e_xy, e_xz], [r_xy, r_xz]], info


def animate(arrays, info):
    # Animate e and r as two subplots of imshow (image representation of a 2D mesh)
    fig, ax = plt.subplots(2, 2)
    header = f"3D Mesh={int(info['h'])}*{int(info['w'])}*{int(info['d'])} elements, dx={1000*info['dx']} mm\n"
    fig.suptitle(header + f"time={int(1000*info['times'][0]):3d} ms", fontsize=12)
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
        fig.suptitle(header + f"time={int(1000*info['times'][num]):3d} ms", fontsize=12)
        for i in range(2):
            for j in range(2):
                ax[i, j].imshow(arrays[i][j][num], vmin=0.0, vmax=1.0)

    ani = animation.FuncAnimation(fig, new_frame, frames)
    ani.save(f"./result.gif",  fps=10,  writer="imagemagick")


if __name__ == "__main__":
    arrays, info = read_files()
    animate(arrays, info)
