import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.gridspec import GridSpec
from matplotlib.legend_handler import HandlerPatch
from matplotlib.ticker import MaxNLocator, MultipleLocator
import matplotlib.lines as mlines
import seaborn as sns
import argparse
import warnings
import math


def seconds_to_hms(seconds):
    seconds %= 86400
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    return f"{int(hours):02d}:{int(minutes):02d}"


def hmsxtick(get_long, lines, maxepnum):
    pattern_tot = r'=TOTSUM\s+(\d{4}/\d{2}/\d{2})\s+(\d{2}:\d{2}:\d{2})\s+(\d{4}/\d{2}/\d{2})\s+(\d{2}:\d{2}:\d{2})\s+\d+\.\d+\s+(\d+)'
    for j in range(get_long):
        line = lines[j]
        matches_tot = re.findall(pattern_tot, line)
        if matches_tot:
            match = matches_tot[0]
            begin_date = match[0]
            begin_time = match[1]
            end_date = match[2]
            end_time = match[3]
            sample = int(match[4])
            split_result = begin_time.split(":")
            begin_h = int(split_result[0])
            begin_m = int(split_result[1])
            begin_s = int(float(split_result[2]))  # Converts seconds part to floating point and then to integer
            split_result = end_time.split(":")
            end_h = int(split_result[0])
            end_m = int(split_result[1])
            end_s = int(float(split_result[2]))  # Converts seconds part to floating point and then to integer
            break
    epnum_xtick = np.linspace(0, maxepnum, 7, dtype=int) * sample + (begin_h * 3600 + begin_m * 60 + begin_s)
    hms_xtick = [seconds_to_hms(seconds) for seconds in epnum_xtick]
    return hms_xtick


def ecef_to_geodetic(ecef_coords):
    """
    Conversion of ECEF coordinates to geographical coordinates (latitude, longitude and altitude)

    Args:
        ecef_coords (list or tuple): list or tuple of length 3 containing ECEF coordinates [x, y, z] in meters

    Returns:
        list: list containing longitudes, latitudes, and altitudes (lon, lat, h)
    """

    # Checks if the input is a valid list or tuple of length 3.
    if not isinstance(ecef_coords, (list, tuple)) or len(ecef_coords) != 3:
        raise ValueError("Input must be a list or tuple of length 3")
    x, y, z = ecef_coords

    # Check if the coordinate value is infinity or NaN
    if any(math.isinf(coord) or math.isnan(coord) for coord in ecef_coords):
        raise ValueError("Input coordinates must be finite numbers")

    # Check if [0, 0, z]
    if x == 0 and y == 0:
        raise ValueError("The coordinates are [0,0,z]")

    # WGS84 Ellipsoid Parameters
    a = 6378137.0  # Long half shafts (meters)
    f = 1 / 298.257223563  # flatness
    b = a * (1 - f)  # bail (arched pot handle)
    e_sq = f * (2 - f)  # First eccentricity squared
    e_prime_sq = e_sq / (1 - e_sq)  # Second eccentricity squared

    x, y, z = ecef_coords

    # Calculating Longitude
    lon = math.atan2(y, x)

    # Calculate the parameter p
    p = math.sqrt(x ** 2 + y ** 2)

    # Initial guess of latitude
    lat = math.atan2(z, p * (1 - e_sq))
    lat_prev = 0
    tolerance = 1e-12  # convergence tolerance
    max_iterations = 100
    iteration = 0

    while abs(lat - lat_prev) > tolerance and iteration < max_iterations:
        lat_prev = lat
        sin_lat = math.sin(lat)
        N = a / math.sqrt(1 - e_sq * sin_lat ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e_sq * (N / (N + h))))
        iteration += 1

    # Calculated Height
    sin_lat = math.sin(lat)
    N = a / math.sqrt(1 - e_sq * sin_lat ** 2)
    h = p / math.cos(lat) - N

    # Converting Longitude and Latitude to Degrees
    lon_deg = math.degrees(lon)
    lat_deg = math.degrees(lat)

    return [lon_deg, lat_deg, h]


def drawskyel(input_file, output_path):
    """
        Plots the satellite sky map and elevation angle distribution based on the input QC file and saves two images.

        This function reads a QC file to generates two images:
        a sky map and an elevation angle distribution map.
        The generated images will be saved with names derived from the QC file.
        For example, if the input QC file is 'x.qc', the resulting images will be saved as 'x_skymap.png' and 'x_elmap.png'.

        Args:
            input_file (str): The path to the QC file containing satellite data.
            output_path (str): The directory where the generated images will be saved.

        Returns:
            str: Returns 'ok' if the images are successfully generated. If an error occurs,
            it returns a string describing the error.

        Raises:
            FileNotFoundError: If the QC file does not exist.
            ValueError: If the QC file is in an incorrect format.

        Example:
            # >>> drawskyel("path/to/file.qc", "output/directory/")
            # This will generate two images: 'file_skymap.png' and 'file_elmap.png' in the 'output/directory/'.
    """

    try:
        with open(input_file, 'r') as file:
            fileContents = file.read()
    except FileNotFoundError:
        return f"File '{input_file}' not found."
    except Exception as e:
        return f"An error occurred: {e}"

    # read filename as sitename
    sitename_match = re.search(r'(.+?)\.qc', os.path.basename(input_file))
    if sitename_match:
        sitename = sitename_match.group(1)
    else:
        return "Not a QC file"
    # create output path of png
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Split file contents into lines
    lines = fileContents.split('\n')
    # Initializing Storage Variables
    # data = []
    sys_data = {
        "G": [],
        "R": [],
        "E": [],
        "C": [],
        "J": []
    }
    signal = {
        "GPS": "",
        "GLO": "",
        "GAL": "",
        "BDS": "",
        "QZS": ""
    }
    # Defining the system
    sys = ['G', 'R', 'E', 'C', 'J']
    sysname = ['GPS', 'GLO', 'GAL', 'BDS', 'QZS']
    file_sys = []
    # Use regular expressions to get L1-signal
    pattern_begin = r'#SYS Freq.Band'
    pattern_end = r'#Sky  prn'
    begin = None
    end = None
    for j, line in enumerate(lines):
        if re.search(pattern_begin, line):
            begin = j + 1
            break
    for j, line in enumerate(lines):
        if re.search(pattern_end, line):
            end = j - 1
            break
    if begin == None:
        return 'no sig data'
    if end == None:
        return 'no sky data'
    for j in range(begin, end):
        line = lines[j]
        parts = line.split()
        sig_sys = parts[0][1:]
        if signal[sig_sys] == "":
            signal[sig_sys] = f"S{parts[2]}"
    # Use regular expressions to match and extract the required data
    pattern_arc = r'arc\s+([a-zA-Z0-9]+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)'
    pattern_normal = r'\s+([a-zA-Z0-9]+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)'
    for line in lines:
        matches_arc = re.findall(pattern_arc, line)
        matches_normal = re.findall(pattern_normal, line)
        if matches_arc:
            match = matches_arc[0]
            dataRow = (
                1,  # mark arc
                match[0],
                int(match[1]),
                float(match[2]),
                float(match[3]),
                int(match[4])
            )
            if match[0][0] in sys and dataRow[-1] > 0:
                sys_data[match[0][0]].append([float(match[2]), float(match[3]), int(match[4]), int(match[1])])
                # data.append(dataRow)
                if match[0][0] not in file_sys:
                    file_sys.append(match[0][0])
        elif matches_normal:
            match = matches_normal[0]
            dataRow = (
                0,  # mark normal
                match[0],
                int(match[1]),
                float(match[2]),
                float(match[3]),
                int(match[4])
            )
            if match[0][0] in sys and dataRow[-1] > 0:
                sys_data[match[0][0]].append([float(match[2]), float(match[3]), int(match[4]), int(match[1])])
                # data.append(dataRow)
                if match[0][0] not in file_sys:
                    file_sys.append(match[0][0])
    if len(file_sys) == 0:
        return 'no sky data'
    # sky = pd.DataFrame(data, columns=['arc', 'prn', 'epnum', 'az', 'el', 'snr'])

    # Number of systems
    sysnum = len(file_sys)
    col = (sysnum + 1) // 2

    # plot skymap
    for k in range(sysnum):
        az = [np.radians(item[0]) for item in sys_data[file_sys[k]]]
        el = [90 - item[1] for item in sys_data[file_sys[k]]]
        snr = [item[2] for item in sys_data[file_sys[k]]]
        current_sys = sysname[sys.index(file_sys[k])]

        ax = plt.subplot(2, col, k + 1, projection='polar')
        ax.plot(np.radians([0]), [0], 'rp', markersize=10)
        ax.set_rmax(90)
        ax.set_rticks([0, 15, 30, 45, 60, 75, 90])
        ax.tick_params(axis='x', pad=10)
        ax.set_yticklabels(['', '75', '60', '45', '30', '15', ''], fontsize=16)
        ax.set_thetagrids(np.arange(0, 360, 30), fontsize=16)
        ax.set_theta_direction(-1)  # Set the angle to increment counterclockwise
        ax.set_theta_offset(np.radians(90))  # Place the 0° at the top
        plt.title(f'{current_sys} SkyMap', fontsize=18)

        scatter = ax.scatter(az, el, c=snr, cmap='plasma', vmin=min(snr), vmax=max(snr), alpha=0.75, s=10)
        # add color bar to show the correspondence between SNR and color
        cbar = plt.colorbar(scatter, ax=ax, pad=0.12)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label(f'{current_sys} {signal[current_sys]}', fontsize=16)
        # Automatic setting of integer scales
        cbar.locator = MaxNLocator(integer=True)  # Make the color bar scale an integer
        cbar.update_ticks()  # Updating the scale to apply settings

    # Output settings
    plt.suptitle(f'{sitename} SkyMap', fontsize=24, ha='center', y=0.91)
    plt.subplots_adjust(top=0.85, wspace=0.15)
    output_file = os.path.join(output_path, f'{sitename}_skymap.png')
    # If the output file already exists, delete the original file
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.gcf().set_size_inches(4 * sysnum, 14)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    # plot elmap
    for k in range(sysnum):
        epnum = [item[-1] for item in sys_data[file_sys[k]]]
        el = [90 - item[1] for item in sys_data[file_sys[k]]]
        snr = [item[2] for item in sys_data[file_sys[k]]]
        current_sys = sysname[sys.index(file_sys[k])]

        plt.subplot(2, col, k + 1)
        plt.ylim(0, 90)
        plt.yticks(np.arange(0, 91, 10), fontsize=14)
        plt.xlim(0, max(epnum))
        hms_xtick = hmsxtick(10, lines, max(epnum))
        plt.xticks(np.linspace(0, max(epnum), 7), hms_xtick, fontsize=14)
        plt.xlabel('Time', fontsize=16)
        plt.ylabel('Elevation Angle/°', fontsize=16)
        plt.tick_params(axis='x', direction='in', pad=8)
        plt.tick_params(axis='y', direction='in', pad=6)
        plt.title(f'{current_sys} ElMap', fontsize=18)

        scatter = plt.scatter(epnum, el, c=snr, cmap='plasma', vmin=min(snr), vmax=max(snr), alpha=0.75, s=10)
        # add color bar
        cbar = plt.colorbar(scatter)
        cbar.set_label(f'{current_sys} {signal[current_sys]}', fontsize=14)
        cbar.ax.tick_params(labelsize=14)

    plt.suptitle(f'{sitename} ElMap', fontsize=24, ha='center', y=0.9)
    plt.subplots_adjust(top=0.85, wspace=0.088)
    output_file = os.path.join(output_path, f'{sitename}_elmap.png')
    # If the output file already exists, delete the original file
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.gcf().set_size_inches(4 * sysnum, 14)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    return 'ok'


def drawpos(input_file, output_path):
    """
        Plots the distribution of DOPs (GDOP, PDOP, HDOP, VDOP) and ENU time-series of SPP
        from a QC file and saves the images as x_dop.png and x_enu.png.

        This function reads a QC file that contains DOP data (GDOP, PDOP, HDOP, VDOP) and SPP ENU results and
        generates plots showing the distribution of these values and time-series of SPP ENU. The generated images
        will be saved as 'x_dop.png' and 'x_enu.png', where 'x' is derived from the input QC file name.

        Args:
            input_file (str): The path to the QC file containing DOP data.
            output_path (str): The directory where the generated image will be saved.

        Returns:
            str: Returns 'ok' if the image is successfully generated. If an error occurs, it returns a string
            describing the error.

        Raises:
            FileNotFoundError: If the QC file does not exist.
            ValueError: If the QC file is in an incorrect format or missing DOP data.

        Example:
            # >>> result = drawpos("path/to/file.qc", "output/directory/")
            # This will generate images 'file_dop.png' and 'file_enu.png' in the 'output/directory/'.
    """
    try:
        with open(input_file, 'r') as file:
            fileContents = file.read()
    except FileNotFoundError:
        return f"File '{input_file}' not found."
    except Exception as e:
        return f"An error occurred: {e}"

    # read filename as sitename
    sitename_match = re.search(r'(.+?)\.qc', os.path.basename(input_file))
    if sitename_match:
        sitename = sitename_match.group(1)
    else:
        return "Not a QC file"
    # create output dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file = os.path.join(output_path, f'{sitename}_dops.png')
    # Split file contents into lines
    lines = fileContents.split('\n')
    # Initializing Storage Variables
    data = []
    ave_spp = [0.0, 0.0, 0.0]
    # Use regular expressions to match and extract the required data
    pattern_begin = r'#POS'
    pattern_end = r'#list'
    begin = None
    finish = None
    for j, line in enumerate(lines):
        if re.search(pattern_begin, line):
            begin = j
            coords = re.findall(r'[-]?\d+\.\d+', line)
            ave_spp[0], ave_spp[1], ave_spp[2] = map(float, coords)
        if re.search(pattern_end, line):
            finish = j
            break
    if begin == None:
        return 'no pdop data'
    if finish == None:
        finish = len(lines)
    pattern = r'(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)'
    for j in range(begin, finish):
        match = re.search(pattern, lines[j])
        if match:
            enu = [float(match.group(6)), float(match.group(7)), float(match.group(8))]
            dataRow = (int(match.group(1)), float(match.group(2)), float(match.group(3)),
                       float(match.group(4)), float(match.group(5)), enu)
            data.append(dataRow)
    if len(data) == 0:
        return 'no dop data or format error'
    # save data
    POS = pd.DataFrame(data, columns=['epnum', 'gdop', 'pdop', 'hdop', 'vdop', 'spp_enu'])

    # Plot distribution of DOPs
    # Setting the vertical coordinate range
    min_dop = np.min([np.min(POS['gdop']), np.min(POS['pdop']), np.min(POS['hdop']), np.min(POS['vdop'])])
    max_dop = np.max([np.max(POS['gdop']), np.max(POS['pdop']), np.max(POS['hdop']), np.max(POS['vdop'])])
    if np.floor(min_dop / 0.1) == 0:
        y1 = 0
    else:
        y1 = (np.floor(min_dop / 0.1) - 1) * 0.1
    y2 = np.ceil(max_dop / 0.1) * 0.1
    plt.ylim(y1, y2)
    yticks = np.linspace(y1, y2, 10)
    plt.yticks(yticks, labels=[f"{val:.2f}" for val in yticks], fontsize=14)
    # Set the x-axis to time
    plt.xlim(0, np.max(POS['epnum']))
    hms_xtick = hmsxtick(10, lines, np.max(POS['epnum']))
    plt.xticks(np.linspace(0, np.max(POS['epnum']), 7), hms_xtick, fontsize=14)
    # Setting the label
    plt.xlabel('Time', fontsize=16)
    plt.ylabel('DOPs', fontsize=16)
    # Set the x-axis scale to face inward
    plt.tick_params(axis='x', direction='in', pad=8)
    # Set the y-axis scale to face inward
    plt.tick_params(axis='y', direction='in', pad=6)
    # plot
    plt.plot(POS['epnum'], POS['pdop'], linewidth=2, color='#0A6DAF', label='PDOP')
    plt.plot(POS['epnum'], POS['gdop'], linewidth=2, color='#A72061', label='GDOP', linestyle='--')
    plt.plot(POS['epnum'], POS['hdop'], linewidth=2, color='#926AAD', label='HDOP', linestyle='-.')
    plt.plot(POS['epnum'], POS['vdop'], linewidth=2, color='#059244', label='VDOP', linestyle=':')

    # Add grid with customization
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.legend(loc='best', fontsize=12)
    plt.title(f'{sitename} DOPs distribution', fontsize=20, loc='center')
    plt.gcf().set_size_inches(16, 10)
    # If the output file already exists, delete the original file
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    # Plot SPP ENU time-series
    # Create a 3x1 subplot with a shared x-axis
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 6))
    # Label and color settings
    labels = ['E-W (m)', 'N-S (m)', 'U-D (m)']
    color_list = ['#0070b8', '#fbb01a', '#d40d8c']
    epnum = POS['epnum']

    for i in range(3):
        enu = [lst[i] for lst in POS['spp_enu']]
        # 计Calculate mean, standard deviation and root mean square error
        ave = np.mean(enu)
        std = np.std(enu)
        rms = np.sqrt(np.mean(np.square(enu)))
        # Plotting each subgraph
        axes[i].plot(epnum, enu, '.-', linewidth=2, color=color_list[i])
        axes[i].set_ylabel(labels[i], fontsize=16, labelpad=5)
        # Setting the scale size
        axes[i].tick_params(axis='both', labelsize=14, direction='in', pad=6)
        # Add a horizontal line to the 0-axis
        axes[i].axhline(0, color='black', linewidth=1)
        # Setting the gridlines
        axes[i].grid(True, linestyle='--')
        # Shows AVE, STD and RMS in the upper right corner
        textstr = f'AVE={ave:.4f}m STD={std:.4f}m RMS={rms:.4f}m'
        text_height = 0.9 if i == 0 else 0.97
        axes[i].text(0.99, text_height, textstr, transform=axes[i].transAxes,
                     fontsize=12, verticalalignment='top', horizontalalignment='right')

    # Only the bottom subplot shows the x-axis scale, the top two hide it
    for ax in axes[:-1]:
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    # Setting the x-axis range and scale labels
    axes[-1].set_xlim(0, max(epnum))
    axes[-1].set_xticks(np.linspace(0, max(epnum), 7))
    axes[-1].set_xticklabels(hms_xtick)
    axes[-1].set_xlabel('Time', fontsize=16, labelpad=4)
    # Add BLH information in the upper right corner
    blh = ecef_to_geodetic(ave_spp)
    lon_str = 'E' if blh[0] >= 0 else 'W'
    lat_str = 'N' if blh[1] >= 0 else 'S'
    axes[0].text(0.99, 0.98, f'ORI= {blh[0]:.9f}°{lon_str} {blh[1]:.9f}°{lat_str} {blh[2]:.4f}m',
                 transform=axes[0].transAxes, fontsize=12,
                 verticalalignment='top', horizontalalignment='right')

    # Setting the title and size
    axes[0].set_title(f'{sitename} SPP-ENU coordinate series', fontsize=20, loc='center', pad=10)
    fig.set_size_inches(16, 10)
    # Adjusting subgraph spacing
    plt.tight_layout()
    # If the output file already exists, delete the original file
    output_file = os.path.join(output_path, f'{sitename}_pos.png')
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    return 'ok'


def drawsys(input_file, output_path):
    """
       Plots the summary of systems from a QC file and saves the image as x_sys.png.

       This function reads a QC file that contains systems qc results and generates a plot showing
       the summary of systems. The generated image will be saved as 'x_sys.png', where 'x' is derived
       from the input QC file name.

       Args:
           input_file (str): The path to the QC file containing DOP data.
           output_path (str): The directory where the generated image will be saved.

       Returns:
           str: Returns 'ok' if the image is successfully generated. If an error occurs, it returns a string
           describing the error.

       Raises:
           FileNotFoundError: If the QC file does not exist.
           ValueError: If the QC file is in an incorrect format or missing DOP data.

       Example:
           # >>> result = drawsys("path/to/file.qc", "output/directory/")
           # This will generate image 'file_sys.png' in the 'output/directory/'.
    """
    try:
        with open(input_file, 'r') as file:
            fileContents = file.read()
    except FileNotFoundError:
        return f"File '{input_file}' not found."
    except Exception as e:
        return f"An error occurred: {e}"

    # read filename as sitename
    sitename_match = re.search(r'(.+?)\.qc', os.path.basename(input_file))
    if sitename_match:
        sitename = sitename_match.group(1)
    else:
        return "Not a QC file"

    # Creating an output folder
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file = os.path.join(output_path, f'{sitename}_sys.png')
    # Split file contents into lines
    lines = fileContents.split('\n')

    # Initializing Storage Variables
    MP = {
        "mp1": [0.0] * 5,
        "mp2": [0.0] * 5,
        "mp3": [0.0] * 5,
        "mp4": [0.0] * 5,
        "mp5": [0.0] * 5,
        "mp6": [0.0] * 5,
        "mp7": [0.0] * 5,
        "mp8": [0.0] * 5,
    }
    data = []
    sys_order = ["GPS", "GLO", "GAL", "BDS", "QZS"]
    # Positioning system information location
    pattern_begin = r'#SYSSUM'
    begin = None
    for j, line in enumerate(lines):
        if re.search(pattern_begin, line):
            begin = j
            finish = min(j + 6, len(lines))
            break
    if begin == None:
        return 'no sys data'

    # Use regular expressions to match and extract the required data
    pattern = r'=([A-Z]+)SUM\s+\d+\s+\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+(\d+)'
    pattern_mp = r'(mp\d+):(\d+\.\d+)'
    for j in range(begin, finish):
        line = lines[j]
        matches = re.findall(pattern, line)
        matches_mp = re.findall(pattern_mp, line)
        if matches:
            match = matches[0]
            if matches_mp:
                for match_mp in matches_mp:
                    MP[match_mp[0]][sys_order.index(match[0])] = float(
                        match_mp[1])  # Storing mp data in a predefined system order
            dataRow = (match[0], float(match[1]), float(match[2]), int(match[3]), int(match[4]))
            data.append(dataRow)
    if len(data) == 0:
        return 'no sys data'
    gnssum = pd.DataFrame(data, columns=['sys', 'rt', 'meansat', 'sig', 'oslps'])
    sysnum = len(gnssum)
    sysname = gnssum['sys'].tolist()
    indices = [sys_order.index(name) for name in sysname]
    index = np.arange(sysnum)

    # Histogram of completeness
    fig, ax1 = plt.subplots(figsize=(10, 6))
    bars = ax1.bar(index, gnssum['rt'], alpha=0.75, color='#a1d8e8', label='obs_rt', align='center')
    # Plotting the mp stack
    mp_keys = [key for key in MP.keys() if sum(MP[key]) != 0]
    mp_values = [[MP[key][i] for i in indices] for key in mp_keys]
    mpx_sum_max = max(sum(values) for values in mp_values)
    if mpx_sum_max * 100 > min(gnssum['rt']):
        scaler_factor = np.floor(min(gnssum['rt']) / mpx_sum_max)
        scaler = np.floor(scaler_factor / 10) * 10
    else:
        scaler = 100
    mp_values = [[MP[key][i] * scaler for i in indices] for key in mp_keys]
    colors = ['#c85e62', '#f47254', '#f59c7c', '#fcc1a6', '#fded95', '#d0e2c0', '#a1d8e8', '#49c2d9']
    color_list = colors[:len(mp_keys)]
    ax1.stackplot(index, *mp_values, labels=mp_keys, colors=color_list)
    # Plotting the measuresat line graph
    ax1.scatter(index, gnssum['meansat'], marker='*', s=96, color='yellow', label='meansat', linewidth=2)
    # Setting up column axes
    ax1.set_xlabel('Satellite System', fontsize=16)
    mp_unit = "" if scaler == 100 else f"{100 / scaler:.2f}"
    ax1.set_ylabel(f'Obs Completeness (%) / Multipath ({mp_unit}cm)', fontsize=16)
    ax1.set_title(f'{sitename} Systems Summary', fontsize=20, pad=10)
    ax1.set_xticks(index)
    ax1.set_xticklabels(sysname)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.set_ylim(0, 110)
    ax1.grid(True, which='both', axis='y', linestyle='--', linewidth=0.7)
    # Plotting o/slps line graphs
    ax2 = ax1.twinx()
    ax2.plot(index, gnssum['oslps'], marker='o', markersize=8, color='#a2c986', label='o/slip', linewidth=2)
    ax2.set_ylabel('Obs/Cycle-Slip', fontsize=16)
    ax2.set_ylim(0, np.ceil(max(gnssum['oslps']) / 100) * 100 + np.floor(0.08 * max(gnssum['oslps']) / 100) * 100)
    ax2.tick_params(axis='y', labelsize=14)
    ax2.grid(True, which='both', axis='y', linestyle='--', linewidth=0.7)
    # Adding labels on the bars and line
    for i, v in enumerate(gnssum['rt']):
        ax1.text(i, v + 1, f'{v}%', color='black', ha='center', fontsize=12)
    for i, v in enumerate(gnssum['oslps']):
        ax2.text(i - 0.02, v - (max(0.035 * max(gnssum['oslps']), 0)), f'{v}', color='#a2c986', ha='right', fontsize=12)
    for i, v in enumerate(gnssum['meansat']):
        ax1.text(i + 0.03, v - (max(0.035 * max(gnssum['meansat']), 0)), f'{v}', color='yellow', ha='left', fontsize=12)
    # for i in index:
    #     stack_height=0
    #     for mp in mp_values:
    #         if mp[i]==0:
    #             ax1.text(index[i],stack_height+0.1,'0',ha='center', va='center', color='gray')
    #         else:
    #             ax1.text(index[i], stack_height + mp[i] / 2, f'{mp[i]:.1f}', ha='center', va='center')
    #         stack_height += mp[i]
    # Legend and layout adjustments
    fig.tight_layout()
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = labels1 + labels2
    ax1.legend(handles, labels, loc='upper center', ncol=11, fontsize=12, handletextpad=0.3)
    # Output set
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.gcf().set_size_inches(4 * sysnum, 10)
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    return 'ok'


def drawts(input_file, output_file_path):
    # Read the contents of the file
    try:
        with open(input_file, 'r') as file:
            fileContents = file.read()
    except FileNotFoundError:
        return f"File '{input_file}' not found."
    except Exception as e:
        return f"An error occurred: {e}"

    # Retrieve site name
    sitename_match = re.search(r'(.+?)\.qc', os.path.basename(input_file))
    if sitename_match:
        sitename = sitename_match.group(1)
    else:
        return "Not a QC file"

    # Creating an output folder
    if not os.path.exists(output_file_path):
        os.makedirs(output_file_path)
    output_file = os.path.join(output_file_path, f'{sitename}_ts.png')
    # Split file contents into lines
    lines = fileContents.split('\n')

    # Initializing Storage Variables
    data = []

    # Use regular expressions to match and extract the required data
    pattern_begin = r'#list'
    begin = None
    for j, line in enumerate(lines):
        if re.search(pattern_begin, line):
            begin = j
            break
    if begin == None:
        return 'no flag list'
    pattern_arc = r'arc\s+([a-zA-Z0-9]+)\s+(\d+)\s+(\d+)\s+(\d+)'
    pattern = r'\s+([a-zA-Z0-9]+)\s+(\d+)\s+(\d+)\s+(\d+)'

    for j in range(begin, len(lines)):
        line = lines[j]
        matches_arc = re.findall(pattern_arc, line)
        matches = re.findall(pattern, line)
        if matches_arc:
            match = matches_arc[0]
            dataRow = [1, match[0], int(match[1]), int(match[2]), int(match[3])]
            data.append(dataRow)
        if matches:
            match = matches[0]
            dataRow = [0, match[0], int(match[1]), int(match[2]), int(match[3])]
            data.append(dataRow)
    if len(data) == 0:
        return 'no flag list'
    # Defining column names and data types
    column_names = ['arc', 'prn', 'epnum', 'obsflag', 'slipflag']
    data_types = [int, str, int, int, int]  # Specify the data type of each column in column order
    # use Pandas DataFrame to create
    ts = pd.DataFrame(data, columns=column_names)
    ts = ts.astype(dict(zip(column_names, data_types)))  # Converts to a specified data type

    # Stores the starting point of each arc
    index = []
    for j in range(len(ts)):
        if ts['arc'][j] == 1:
            index.append(j)
    index.append(len(ts))

    class FlagData:
        def __init__(self, prn, data):
            self.prn = prn
            self.data = pd.DataFrame(data, columns=['epnum', 'obsflag', 'slipflag'])

    # Storing data by prn
    flag_data = []
    for j in range(len(index) - 1):
        flag1 = FlagData(prn=ts['prn'][index[j]],
                         data=(ts[['epnum', 'obsflag', 'slipflag']][index[j]:index[j + 1]].reset_index(drop=True)))
        flag_data.append(flag1)
    maxepnum = np.max(ts['epnum'])
    hms_xtick = hmsxtick(begin, lines, maxepnum)
    sys = ['G', 'R', 'E', 'C', 'J']

    sys_color_rgb = [[112, 180, 182], [226, 188, 116], [203, 214, 150], [180, 139, 171], [147, 217, 212]]
    sys_color = [[r / 255, g / 255, b / 255] for r, g, b in sys_color_rgb]
    n = 0
    while n < len(flag_data):
        for j in range(flag_data[n].data.shape[0]):
            currentcolor = sys_color[sys.index(flag_data[n].prn[0])]
            # Start - End: current is 1 and next is 2 (j is not the last to ensure there is a next), connecting lines
            if j < flag_data[n].data.shape[0] - 1 and flag_data[n].data['obsflag'][j] == 1 and \
                    flag_data[n].data['obsflag'][j + 1] == 2:
                plt.plot([flag_data[n].data['epnum'][j], flag_data[n].data['epnum'][j + 1]],
                         [n + 1, n + 1], color=currentcolor, linewidth=2.5)
                j += 1
            # There are weekly jumps in the middle section: j is not the last one, the current one is 1 and the next one is 0
            if j < flag_data[n].data.shape[0] - 1 and flag_data[n].data['obsflag'][j] == 1 and \
                    flag_data[n].data['obsflag'][j + 1] == 0:
                k = j + 2
                while k < flag_data[n].data.shape[0]:
                    # Find the first 2 backward. Connect the dots.
                    if flag_data[n].data['obsflag'][k] == 2:
                        plt.plot([flag_data[n].data['epnum'][j], flag_data[n].data['epnum'][k]],
                                 [n + 1, n + 1], color=currentcolor, linewidth=2.5)
                        j = k
                        break
                    k += 1
            # Lone point: j is the last and 1, tracing point
            if j == flag_data[n].data.shape[0] - 1 and flag_data[n].data['obsflag'][j] == 1:
                plt.scatter(flag_data[n].data['epnum'][j], n + 1, s=8, edgecolors=currentcolor, facecolors=currentcolor,
                            marker='o')
            # Lone point: the current j is 1 and the next one is also 1 , tracing point
            elif flag_data[n].data['obsflag'][j] == 1 and flag_data[n].data['obsflag'][j + 1] == 1 and \
                    flag_data[n].data['epnum'][j] != flag_data[n].data['epnum'][j]:
                plt.scatter(flag_data[n].data['epnum'][j], n + 1, s=8, edgecolors=currentcolor, facecolors=currentcolor,
                            marker='o')
        n += 1

    # settings
    prn_ytick = [''] + [flag_data[j].prn for j in range(len(flag_data))]
    plt.xlim([0, maxepnum])
    plt.ylim([0, len(flag_data) + 1])
    plt.xticks(np.linspace(0, maxepnum, 7), hms_xtick, fontsize=14)
    plt.grid(True, linestyle='--')
    plt.gca().yaxis.set_major_locator(MultipleLocator(1))
    yticks = range(len(prn_ytick))  # Generate y-axis scale based on the number of prn_yticks
    plt.gca().set_yticks(yticks)  # Set the y-axis scale first
    plt.gca().set_yticklabels(prn_ytick, fontsize=10)
    plt.xlabel('Time', fontsize=16)
    plt.ylabel('PRN', fontsize=16)
    plt.tick_params(axis='x', direction='in', pad=8)
    plt.tick_params(axis='y', direction='in', pad=6)
    plt.title(f'{sitename} Time series', fontsize=24, ha='center', pad=10)
    plt.gcf().set_size_inches(16, 24)

    # Output
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    return 'ok'


def drawsig(input_file, output_path):
    """
        Plots the results of signals (completeness, mp, noise, snr)
        from a QC file and saves the images as x_sig.png.

        This function reads a QC file that contains the results of signals (completeness, mp, noise, snr)  and
        generates plots showing the distribution of these values. The generated images will be saved as
        'x_sig.png', where 'x' is derived from the input QC file name.

        Args:
            input_file (str): The path to the QC file containing DOP data.
            output_path (str): The directory where the generated image will be saved.

        Returns:
            str: Returns 'ok' if the image is successfully generated. If an error occurs, it returns a string
            describing the error.

        Raises:
            FileNotFoundError: If the QC file does not exist.
            ValueError: If the QC file is in an incorrect format or missing DOP data.

        Example:
            # >>> result = drawpos("path/to/file.qc", "output/directory/")
            # This will generate the image 'file_sig.png'.
    """
    # ===== Part 1: Open file and prepare data for plot =====
    try:
        with open(input_file, 'r') as file:
            fileContents = file.read()
    except FileNotFoundError:
        return f"File '{input_file}' not found."
    except Exception as e:
        return f"An error occurred: {e}"

    # read filename as sitename
    sitename_match = re.search(r'(.+?)\.qc', os.path.basename(input_file))
    if sitename_match:
        sitename = sitename_match.group(1)
    else:
        return "Not a QC file"
    # create output dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_file = os.path.join(output_path, f'{sitename}_sig.png')
    # Split file contents into lines
    lines = fileContents.split('\n')
    # Initializing Storage Variables
    data = []

    # Use regular expressions to match and extract the required data
    pattern_begin = r'#SYS Freq.Band'
    pattern_end = r'#Sky  prn'
    begin = None
    finish = None
    for j, line in enumerate(lines):
        if re.search(pattern_begin, line):
            begin = j
        if re.search(pattern_end, line):
            finish = j
            break
    if begin == None:
        return 'no signal data'
    if finish == None:
        finish = len(lines)
    pattern = (r"^=(\b[A-Z]{3}\b)\s+\w+\s+(\w+)\s+\d+\s+\d+\s+([\d\.]+)\s+\d+\s+\d+\s+([\d\.]+)\s+([\d\.]+)\s+(["
               r"\d\.]+)\s+(\[.*?\])\s+(\[.*?\])")
    for j in range(begin, finish):
        match = re.search(pattern, lines[j])
        if match:
            pnoise = [float(x) for x in re.findall(r"[\d.]+", match.group(7))]
            lnoise = [float(x) for x in re.findall(r"[\d.]+", match.group(8))]
            datarow = (
                match.group(1), match.group(2), float(match.group(3)),
                float(match.group(4)), float(match.group(5)), float(match.group(6)),
                pnoise[:-1], pnoise[-1], lnoise[:-1], lnoise[-1]
            )
            data.append(datarow)
    if len(data) == 0:
        return 'no signal data or format error'
    # save data
    SIG = pd.DataFrame(data, columns=['sys', 'code', 'rt', 'rt10', 'snr', 'mp', 'pbox', 'pnoise', 'lbox', 'lnoise'])
    SIG['sys_code'] = SIG['sys'] + ' ' + SIG['code']

    # Prepare data for map
    snr = np.array(SIG['snr'].tolist()).reshape(-1, 1)
    pnoise_min = min([min(row['pbox']) for i, row in SIG.iterrows()])
    pnoise_max = max([max(row['pbox']) for i, row in SIG.iterrows()])
    lnoise_min = min([min(row['lbox']) for i, row in SIG.iterrows()])
    lnoise_max = max([max(row['lbox']) for i, row in SIG.iterrows()])
    mp_min = min(SIG['mp'])
    mp_max = max(SIG['mp'])
    x_min = min(pnoise_min, lnoise_min, mp_min)
    x_max = max(lnoise_max, pnoise_max, mp_max)
    # Color mapping for bars
    color_bar = ['#7b95c6', '#a2c986', '#fded95', '#f59c7c', '#fcc1a6']
    systems = ['GPS', 'GLO', 'GAL', 'BDS', 'QZS']
    colors = SIG['sys'].map(
        {'GPS': color_bar[0], 'GLO': color_bar[1], 'GAL': color_bar[2], 'BDS': color_bar[3], 'QZS': color_bar[4]})
    # Calculate the difference between rt10 and rt
    rt_diff = SIG['rt10'] - SIG['rt']

    # Create the figure and GridSpec layout
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(1, 2, width_ratios=[1, 3], wspace=0.13)

    # ===== Part 2: Plot bar of signal completeness =====
    ax1 = fig.add_subplot(gs[0])

    # Create a small offset for each bar to manually adjust its position
    y_positions = np.arange(len(SIG['sys_code'])) + 0.5  # 微调位置
    # Set the y-axis range, removing the top and bottom margins
    ax1.set_ylim(min(y_positions) - 0.5, max(y_positions) + 0.5)
    # Plot the main SIG['rt'] bars in reverse order
    ax1.barh(y_positions, -SIG['rt'], color=colors)
    # If rt10 > rt, add the black bars on top of SIG['rt']
    ax1.barh(y_positions, -np.maximum(rt_diff, 0), left=-SIG['rt'], color='black')
    # If rt10 < rt, show the black bars within SIG['rt'] but aligned inside the bar
    ax1.barh(y_positions, np.minimum(rt_diff, 0), color='black')

    # Invert the y-axis so that the bars are drawn from top to bottom
    ax1.invert_yaxis()

    # Y-axis on the right
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.set_yticks(np.arange(len(SIG['sys_code'])) + 0.5)
    ax1.set_yticklabels(SIG['sys_code'], fontsize=10)
    # X-axis settings
    ax1.tick_params(axis='x', direction='inout')
    ax1.set_xlim(-100, 0)
    ax1.set_xticks(np.arange(-100, 1, 20))
    ax1.set_xticklabels(np.arange(100, -1, -20), fontsize=10)
    ax1.set_xlabel('Completeness of Observations/%', fontsize=12)
    # Adding dashed gridlines along the x-axis
    ax1.grid(True, axis='x', linestyle='--', linewidth=0.5)

    # Create bar legend
    # Creating Multicolor Rectangle Handles
    class MultiColorRectHandler(HandlerPatch):
        def __init__(self, **kwargs):
            HandlerPatch.__init__(self, **kwargs)

        def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
            # Create four small rectangular blocks of different colors
            patches_list = []
            total = SIG['sys'].nunique()
            k = 0
            for i, color in enumerate(color_bar):
                if systems[i] not in SIG['sys'].values: continue
                rect = patches.Rectangle([xdescent + k * (width / total), ydescent], width / total, height,
                                         facecolor=color,
                                         transform=trans)
                patches_list.append(rect)
                k += 1
            return patches_list

    # Create a transparent rectangle object to be used as a placeholder for the legend
    rect_color = patches.Rectangle((0, 0), 1, 1, facecolor='none', edgecolor='none')
    # Create a black fill rectangle
    black_rect = patches.Rectangle((0, 0), 1, 1, facecolor='black', edgecolor='black')

    # Create a legend and add custom multi-colored rectangle handles
    ax1.legend(
        [rect_color, black_rect],
        ['Obs_rt', 'Obs_rt(el>10°)'],
        handler_map={rect_color: MultiColorRectHandler()},
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=2,
        fontsize=10
        # handletextpad=0.2,
        # columnspacing=1
    )

    # ===== Part 2: Plot heatmap by snr with noise distribution and mp =====
    ax2 = fig.add_subplot(gs[1])
    # Synchronized y-axis range
    ax2.set_ylim(ax1.get_ylim())
    # plot heatmap
    heatmap = sns.heatmap(snr, annot=False, cmap="YlGnBu", cbar=True, yticklabels=False, xticklabels=False,
                          ax=ax2)
    # Get the axis of the colorbar
    cbar = heatmap.collections[0].colorbar
    # Get the current colorbar position and move it to the right
    cbar_ax = cbar.ax
    pos = cbar_ax.get_position()
    cbar_ax.set_position([pos.x0 + 0.02, pos.y0, pos.width, pos.height])
    # Adjusting the font size of the colorbar label
    cbar_ax.yaxis.label.set_size(10)

    # color set
    color_other = ['#fbb01a', color_bar[2], '#fbb01a', '#f47254', '#f47254', color_bar[3], '#d40d8c', '#49c2d9',
                   '#ff69b4']

    # Hide bottom xticks (including SNR labels)
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    # Labeled SNR
    snr_ax = ax2.twinx()
    snr_ax.set_ylim(ax2.get_ylim())
    snr_ax.set_ylabel('SNR/dB', fontsize=12)
    snr_ax.set_yticks(np.arange(len(SIG['sys_code'])) + 0.5)
    snr_ax.set_yticklabels(SIG['snr'], fontsize=10)

    # Reset the x-axis range for boxplots, making sure it does not stretch the heatmap
    boxplot_ax = ax2.twiny()

    # Plot the box plots on the new x-axis, adjusting only the boxplot x-axis to fit the pnoise distribution
    n = len(SIG)
    for i, row in SIG.iterrows():
        if any(val != 0.0 for val in row['pbox']):
            boxplot_ax.boxplot(row['pbox'], positions=[i + 0.5], widths=0.4, vert=False, patch_artist=True,
                               boxprops=dict(facecolor=color_other[1], color=color_other[0]),  # Customize color
                               whiskerprops=dict(color=color_other[0], linestyle="--", linewidth=1.5),  # Whisker style
                               capprops=dict(color=color_other[0]),  # Cap style
                               medianprops=dict(color=color_other[2], linewidth=2),
                               showfliers=False)
            # plot the average of prange_noise
            # boxplot_ax.plot(row['pnoise'], i + 0.5, 'o',color=color_other[3])

        if any(val != 0.0 for val in row['lbox']):
            boxplot_ax.boxplot(row['lbox'], positions=[i + 0.5], widths=0.4, vert=False, patch_artist=True,
                               boxprops=dict(facecolor=color_other[5], color=color_other[4]),  # Customize color
                               whiskerprops=dict(color=color_other[4], linestyle="--", linewidth=1.5),  # Whisker style
                               capprops=dict(color=color_other[4]),  # Cap style
                               medianprops=dict(color=color_other[6], linewidth=2),  # Median line style
                               showfliers=False)
            # plot the average of carrier_noise
            # boxplot_ax.plot(row['lnoise'], i + 0.5, 'o',color=color_other[7])
        # plot signal mp
        if row['mp'] != 0.0: boxplot_ax.plot(row['mp'], i + 0.5, '*', color=color_other[8], markersize=12)

    # Set the x-axis
    padding = (x_max - x_min) * 0.04  # Add 4% top and bottom padding
    boxplot_ax.xaxis.tick_bottom()
    boxplot_ax.set_xlim(x_min - padding, x_max + padding)
    boxplot_ax.tick_params(axis='x', direction='inout', labelsize=10)
    boxplot_ax.set_xlabel('Noise of Pseudorange or Carrier/m', fontsize=12)
    boxplot_ax.xaxis.set_label_coords(0.5, -0.05)
    # Set the y-axis
    boxplot_ax.set_yticklabels([])
    # Adding dashed gridlines along the x-axis
    boxplot_ax.grid(True, axis='x', linestyle='--', linewidth=0.5)

    # Creating Box Chart Legend Items
    pbox_legend = mlines.Line2D([], [], color=color_other[0], marker='s', markersize=10, linestyle='-',
                                label='Pseudorange Noise Distribution',
                                markerfacecolor=color_other[1], markeredgewidth=2)

    lbox_legend = mlines.Line2D([], [], color=color_other[4], marker='s', markersize=10, linestyle='-',
                                label='Carrier Noise Distribution',
                                markerfacecolor=color_other[5], markeredgewidth=2)
    # Creating a legend for mean points and mp
    # pnoise_legend = mlines.Line2D([], [], color=color_other[3], marker='o', linestyle='None', markersize=8,
    #                               label='Mean of PR Noise')
    # lnoise_legend = mlines.Line2D([], [], color=color_other[7], marker='o', linestyle='None', markersize=8,
    #                               label='Mean of CP Noise')
    mp_legend = mlines.Line2D([], [], color=color_other[8], marker='*', linestyle='None', markersize=10, label='MP')

    # Legend
    boxplot_ax.legend(
        handles=[pbox_legend, lbox_legend, mp_legend],
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=5,
        fontsize=10,
        handletextpad=0.3,  # Control the spacing between graphics and text
        # columnspacing=0.7  # Controls the horizontal spacing between legend items
    )

    # Title
    fig.suptitle(f'{sitename} signal snr, completeness and noise distribution', fontsize=18, y=0.95)
    # Output
    if os.path.exists(output_file):
        os.remove(output_file)
    # save png
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    return 'ok'


def main():
    warnings.filterwarnings("ignore", category=UserWarning)
    parser = argparse.ArgumentParser(description='Plot qc_file')
    parser.add_argument('-f', '--function', choices=['pos', 'sys', 'sig', 'skyel', 'ts', 'all'],
                        help='Function to plot')
    parser.add_argument('-i', '--input', help='Input qc_file')
    parser.add_argument('-o', '--output', help='Output png path')

    args = parser.parse_args()

    if args.function == 'pos':
        res = drawpos(args.input, args.output)
    elif args.function == 'sys':
        res = drawsys(args.input, args.output)
    elif args.function == 'skyel':
        res = drawskyel(args.input, args.output)
    elif args.function == 'ts':
        res = drawts(args.input, args.output)
    elif args.function == 'sig':
        res = drawsig(args.input, args.output)
    elif args.function == 'all':
        res = [drawpos(args.input, args.output), drawsys(args.input, args.output),
               drawskyel(args.input, args.output), drawts(args.input, args.output),
               drawsig(args.input, args.output)]
    else:
        print("Invalid function specified")

    if args.function != 'all':
        # Individual function case
        if res == 'ok':
            print(f'finish plot')
        else:
            print(f'error:{res}')
    else:
        # Multiple function cases ('all')
        if res == ['ok', 'ok', 'ok', 'ok', 'ok']:
            print(f'finish plot')
        else:
            choice = ['pos', 'sys', 'skyel', 'ts', 'sig']
            for i, r in enumerate(res):
                if r != 'ok':
                    print(f'error: {r} in {choice[i]} plot')
                else:
                    print(f'finish {choice[i]} plot')


if __name__ == '__main__':
    main()
