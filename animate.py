import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.ticker as ticker

def read_files(data_path="./data"):
    # Read files from data path, expected on the form 'e_xy_X.csv' 
    # with 'e_' or 'r_', then 'xy_' or 'xz_', then 'X' (number), then '.csv'.
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
    header = f"3D mesh={int(info['h'])}*{int(info['w'])}*{int(info['d'])} elements, time="
    title_size = 10
    subtitle_size = 8
    emax = 1.0
    rmax = 0.2
    fig.suptitle(header + f"{float(info['times'][0]):1.2f}", fontsize=title_size)
    plt.rc('axes', titlesize=subtitle_size)
    e_xy = arrays[0][0]
    frames = len(e_xy)
    titles = [["e front (xy)", "e right (xz)"], ["r front (xy)", "r right (xz)"]]
    for i in range(2): # e or r
        for j in range(2): # xy or xz
            im = ax[i, j].imshow(arrays[i][j][0], vmin=0.0, vmax=emax if (i == 0) else rmax)
            ax[i, j].set_title(titles[i][j])
            ax[i, j].xaxis.set_major_locator(ticker.NullLocator())
            ax[i, j].yaxis.set_major_locator(ticker.NullLocator())
            fig.colorbar(im, ax=ax[i][j])

    # For matplotlibs FuncAnimation to call for each new frame
    def new_frame(num):
        fig.suptitle(header + f"{float(info['times'][num]):1.2f}", fontsize=title_size)
        for i in range(2): # e or r
            for j in range(2): # xy or xz
                ax[i, j].imshow(arrays[i][j][num], vmin=0.0, vmax=emax if (i == 0) else rmax)

    ani = animation.FuncAnimation(fig, new_frame, frames)
    ani.save(f"./result.gif",  fps=10,  writer="imagemagick")


if __name__ == "__main__":
    arrays, info = read_files()
    animate(arrays, info)
