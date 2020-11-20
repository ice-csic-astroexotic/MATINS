#!/usr/bin/evn python3
# -*- coding: utf-8 -*-

""" Regression test for Magneto-Thermal 2D.

    Authors:

        Alberto Garcia Garcia (garciagarcia@ice.csic.es)

    Copyright (c) MAGNESIA (ICE-CSIC)

"""

import argparse
import json
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import typing
import urllib.error
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

log = logging.getLogger(__name__)

def configure_plotting() -> None:
    """
    Configures Matplotlib font sizes.
    """

    plt.rc('font', size=24) # Controls default text sizes.
    plt.rc('axes', titlesize=32) # Fontsize of the axes title.
    plt.rc('axes', labelsize=32) # Fontsize of the x and y labels.
    plt.rc('xtick', labelsize=32) # Fontsize of the tick labels.
    plt.rc('ytick', labelsize=32) # Fontsize of the tick labels.
    plt.rc('legend', fontsize=24) # Legend fontsize.
    plt.rc('figure', titlesize=32) # Fontsize of the figure title.

def read_profile(
        filename: pathlib.Path
) -> typing.Tuple[np.array, np.array, np.array]:
    """
    Reads a ygraph profile file.

    Args:
        filename: the file name of the profile.

    Returns:
        A tuple of three arrays:
          2D array with one row for each time stamp in the profile.
          1D array with all the timestamps.
          2D array with all coordinates for the profiles for each timestamp.

    """

    profile: list = []
    coordinates: list = []
    timestamps: list = []

    with open(filename, "r") as file:

        current_timestamp: float = 0.0
        current_coordinates: list = []
        current_profile: list = []

        for line in file:

            split = line.strip().split()

            if line == "\n":

                profile.append(np.array(current_profile))
                current_profile = []
                timestamps.append(current_timestamp)
                current_timestamp = 0.0
                coordinates.append(np.array(current_coordinates))
                current_coordinates = []

            elif split[0].split('=')[0] == "\"Time":

                current_timestamp = float(split[1])
                continue

            elif split[0].split('=')[0] == "\"Label":

                continue

            else:

                current_coordinates.append(float(split[0]))
                current_profile.append(float(split[1]))


    return np.array(profile), np.array(timestamps), np.array(coordinates)


def read_magnetic_field(
        filename: pathlib.Path
) -> typing.Tuple[np.array, np.array, np.array, np.array]:
    """
    Read the magnetic field evolution output from the magneto-thermal 2D code.

    Args:
        filename: the file name of the magnetic field map.
        max_snapshots: maximum number of snapshots to read.

    Returns:

    """

    polar_coordinates: list = []
    bphi: list = []
    aphi: list = []
    times: list = []

    # Using a for loop to avoid memory problems when loading a big file.
    with open(filename, "r") as file:

        snapshot = 0

        current_bphi: list = []
        current_aphi: list = []
        current_polar_coordinates: list = []

        for line in file:

            split = line.strip().split()

            if split[0] == "snapshot":

                times.append(float(split[1]))

            elif split[0] == "snapshot_end":

                snapshot += 1

                bphi.append(current_bphi)
                aphi.append(current_aphi)
                polar_coordinates.append(current_polar_coordinates)

                current_bphi = []
                current_aphi = []
                current_polar_coordinates = []

            else:

                current_polar_coordinates.append(
                    [float(split[0]),
                     float(split[1])]
                )
                current_bphi.append(float(split[4]))
                current_aphi.append(
                    float(split[5])
                    * float(split[1])
                    * np.sin(float(split[0]))
                )

    polar_coordinates = np.array(polar_coordinates)
    bphi = np.array(bphi)
    aphi = np.array(aphi)
    times = np.array(times)

    return polar_coordinates, bphi, aphi, times


def read_temperature(
        filename: pathlib.Path,
) -> typing.Tuple[np.array, np.array, np.array]:
    """
    Read the temperature map file from the magneto-thermal 2D evolution code.

    Args:
        filename: the file name of the temperature map.
        max_snapshots: maximum number of snapshots to read.

    Returns:
        A tuple of np.array that contains the polar coordinates of the grid,
        the temperatures associated with each pair of coordinates, and the
        timesteps for each snapshot.
    """

    polar_coordinates: list = []
    temperatures: list = []
    times: list = []

    # Using a for loop to avoid memory problems when loading a big file.
    with open(filename, "r") as file:

        current_polar_coordinates: list = []
        current_temperatures: list = []

        for line in file:

            split = line.strip().split()

            if split[0] == "label":
                continue
            elif split[0] == "snapshot":
                times.append(float(split[1]))
            elif split[0] == "snapshot_end":
                polar_coordinates.append(current_polar_coordinates)
                temperatures.append(current_temperatures)
                current_polar_coordinates = []
                current_temperatures = []
            else:
                current_polar_coordinates.append(
                    [
                        float(split[0]),
                        float(split[1])
                    ]
                )
                current_temperatures.append(float(split[2]))

    polar_coordinates = np.array(polar_coordinates)
    temperatures = np.array(temperatures)
    times = np.array(times)

    return polar_coordinates, temperatures, times


def read_grid_size(
        input_path: pathlib.Path
) -> typing.Tuple[int, int]:
    """
    Grid size reader.

    Extracts the size of the grid (kmax, lmax) from the input file.

    Args:
        input_path: Path for the input file.

    Returns:
        A tuple containing the grid dimensions [kmax, lmax].

    """

    kmax = -1
    lmax = -1

    with open(input_path) as file:

        lines = file.read().splitlines()
        # In the actual file format, the kmax and lmax parameters are located
        # in the third row.
        values = lines[2].split()
        kmax = int(values[0])
        lmax = int(values[1])

    return kmax, lmax

def plot_profile_error(
        gold_profile: np.array,
        run_profile: np.array,
        output_path: pathlib.Path
):

    profile_error = np.abs(gold_profile - run_profile)

    fig = plt.figure(figsize=(32, 16))
    ax = fig.add_subplot(111)
    plot = ax.pcolor(
        profile_error,
        vmin=0.0,
        vmax=np.max(profile_error),
        cmap='inferno',
        edgecolors='k',
        linewidths=1
    )
    ax.set_title(output_path.name)
    fig.colorbar(plot, ax=ax)
    fig.savefig(output_path)
    plt.close()


def plot_profile_errors(
        gold_profiles: np.array,
        run_profiles: np.array,
        coordinates: np.array,
        timestamps: np.array,
        output_path: pathlib.Path
) -> None:

    num_profiles = len(gold_profiles)

    fig = plt.figure(figsize=(32, 16))
    fig.tight_layout()

    for i, profile in enumerate(gold_profiles):

        profile_error = np.abs(gold_profiles[profile] - run_profiles[profile])

        ax = fig.add_subplot(1, num_profiles, i+1)
        plot = ax.pcolor(
            profile_error,
            vmin=0.0,
            vmax=np.max(profile_error),
            cmap='inferno',
        )
        ax.set_title(profile)
        fig.colorbar(plot, ax=ax)

    fig.savefig(output_path)
    plt.close()

def plot_error(
        angles: np.array,
        radial: np.array,
        absolute_error: np.array,
        relative_error: np.array,
        units: str,
        output_filename: pathlib.Path
):
    """
    Error plot generator.

    Generates a heatmap-like plot of both the absolute and relative error
    in two semicircles in polar coordinates.

    Args:
        angles: Set of angular divisions for the plot.
        radial: Set of radial divisions for the plot.
        absolute_error: Absolute error for each polar coordinate (angle-major).
        relative_error: Relative error for each polar coordinate (angle-major).
        units: Units for the absolute error.
        output_filename: Where to generate the plot and save it.

    Returns:
        Nothing, but generates a figure in the specified path.

    """

    angular_divisions = len(angles)
    radial_divisions = len(radial)

    fig = plt.figure(figsize=(32, 16))
    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.12, top=0.85)
    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    spec = matplotlib.gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
    ax1 = fig.add_subplot(spec[0, 0], polar=True)

    ax1.set_title("", pad=32)

    error_min = np.min(absolute_error)
    error_max = np.max(absolute_error) + 0.001
    error_levels = np.linspace(error_min, error_max, 64)

    # Temperature plot.
    contour_absolute_error = ax1.contourf(
        angles - np.pi / 2,
        radial,
        absolute_error.reshape(
            (angular_divisions, radial_divisions)
        ).transpose(),
        levels=error_levels,
        cmap="hot"
    )

    # Setup colorbar for temperature plot.
    axin1 = inset_axes(
        ax1,
        width='5%',
        height='100%',
        loc='right',
        bbox_to_anchor=(0.15, 0, 1, 1),
        bbox_transform=ax1.transAxes
    )
    colorbar_absolute_error = plt.colorbar(
        contour_absolute_error,
        cax=axin1,
        orientation='vertical'
    )
    colorbar_absolute_error.set_label(
        r"Absolute Error $ {} $".format(units), labelpad=64
    )
    colorbar_absolute_error.set_ticks(
        [
            np.linspace(
                error_min,
                error_max,
                num=5,
                endpoint=True
            )
        ]
    )
    colorbar_absolute_error.set_ticklabels(
        [
            "${:0.2f}$".format(x) for x in colorbar_absolute_error.get_ticks()
        ]
    )

    error_relative_min = np.min(relative_error)
    error_relative_max = np.max(relative_error) + 0.0001
    error_relative_levels = np.linspace(
        error_relative_min,
        error_relative_max,
        64
    )

    # Magnetic field plot.
    contour_relative_error = ax1.contourf(
        angles + np.pi / 2,
        radial,
        relative_error.reshape(
            (angular_divisions, radial_divisions)
        ).transpose(),
        levels=error_relative_levels,
        cmap="inferno"
    )

    # Setup colorbar for magnetic field plot.
    axin2 = inset_axes(
        ax1,
        width='5%',
        height='100%',
        loc='center left',
        bbox_to_anchor=(-0.15, 0, 1, 1),
        bbox_transform=ax1.transAxes
    )
    colorbar_relative_error = plt.colorbar(
        contour_relative_error,
        cax=axin2,
        orientation='vertical'
    )
    colorbar_relative_error.set_label(r"Relative Error", labelpad=-256)
    colorbar_relative_error.set_ticks(
        [
            np.linspace(
                error_relative_min,
                error_relative_max,
                num=5,
                endpoint=True
            )
        ]
    )
    colorbar_relative_error.ax.yaxis.set_ticks_position('left')

    fig.savefig(output_filename)


def relative_l2_norm(
        gold_values: np.array,
        run_values: np.array
) -> float:
    """
    Compute relative L2 norm.

    Computes the relative L2 norm of two arrays of values.

    Args:
        gold_values: Reference values to compare.
        run_values: Run values.

    Returns:
        The relative L2 error between the the arrays.

    """

    sum_squared_reference = np.sum(np.square(gold_values))
    squared_differences = np.sum(np.square(gold_values - run_values))
    relative_l2 = squared_differences / sum_squared_reference

    return relative_l2


def read_regression_configuration(
        filename: str
) -> dict:
    """
    Regression test configuration reader.

    Reads a JSON file with the full regression test configuration and generates
    a dictionary containing such configuration.

    Args:
        filename: File name for the JSON configuration file.

    Returns:
        A dictionary containing the configuration for the test suite.

    """

    log.info("Reading regression configuration file {}".format(filename))

    with open(filename) as file:
        configuration = json.load(file)

    return configuration


def checkout_file(
        file_path: pathlib.Path,
        file_url: str
) -> None:
    """
    Checkout a required test file from a remote server.

    Checks if the requested test file exists in the appropriate path in the
    system. If the path does not exist, it creates it. If the files does not
    exist, it tries to download it from the specified url.

    Args:
        file_path: complete file path (parent + filename) of the test file.
        file_url: remote server URL to download the file if needed.

    Returns:
        Nothing.

    """

    log.info("Checking required test file {}".format(file_path))

    if not file_path.parent.exists():

        file_path.parent.mkdir(parents=True, exist_ok=False)

    if not file_path.exists():

        log.info("File {} needed for test not found...".format(file_path))

        try:

            log.info("Trying to download file...")
            urllib.request.urlretrieve(file_url, file_path)
            log.info("Downloaded!")

        except urllib.error.HTTPError as exception:

            log.error("Download failed: {}".format(exception))

        except urllib.error.URLError as exception:

            log.error("Download failed: {}".format(exception))


def checkout_test_files(
        test_configuration: dict
) -> None:
    """
    Checkout test files.

    Check if all the needed files for the current test exist: magnetic field
    map, temperature map, and all execution input files. If they do not,
    try to fetch them online.

    Args:
        test_configuration: Dictionary with all test configuration.

    Returns:
        Nothing.

    """

    # Checkout input files.
    for input_file, input_file_url in zip(
            test_configuration["input_files"],
            test_configuration["input_urls"]
    ):

        checkout_file(
            pathlib.Path().joinpath(
                "gold_runs",
                test_configuration["folder"],
                input_file
            ),
            input_file_url
        )

    # Checkout magnetic field file.
    checkout_file(
        pathlib.Path().joinpath(
            "gold_runs",
            test_configuration["folder"],
            test_configuration["magnetic_field"]
        ),
        test_configuration["magnetic_field_url"]
    )

    # Checkout temperature map file.
    checkout_file(
        pathlib.Path().joinpath(
            "gold_runs",
            test_configuration["folder"],
            test_configuration["temperature_map"]
        ),
        test_configuration["temperature_map_url"]
    )

    # Checkout magnetic radial profiles.
    if "magnetic_radial_profiles" in test_configuration:
        for profile, profile_url in zip(
                test_configuration["magnetic_radial_profiles"],
                test_configuration["magnetic_radial_profiles_urls"]
        ):
            checkout_file(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                ),
                profile_url
            )

    # Checkout magnetic meridional profiles.
    if "magnetic_meridional_profiles" in test_configuration:
        for profile, profile_url in zip(
                test_configuration["magnetic_meridional_profiles"],
                test_configuration["magnetic_meridional_profiles_urls"]
        ):
            checkout_file(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                ),
                profile_url
            )

    # Checkout temperature radial profiles.
    if "temperature_radial_profiles" in test_configuration:
        for profile, profile_url in zip(
                test_configuration["temperature_radial_profiles"],
                test_configuration["temperature_radial_profiles_urls"]
        ):
            checkout_file(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                ),
                profile_url
            )

    # Checkout temperature meridional profiles.
    if "temperature_meridional_profiles" in test_configuration:
        for profile, profile_url in zip(
                test_configuration["temperature_meridional_profiles"],
                test_configuration["temperature_meridional_profiles_urls"]
        ):
            checkout_file(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                ),
                profile_url
            )


def create_test_execution_folders(
        test_configuration: dict
) -> None:
    """
    Creates all test execution folders.

    Creates all the needed folders for the execution of the test: the results,
    runs, and all execution input and output folders (in, out, and outb).

    Args:
        test_configuration: Dictionary with all test configuration.

    Returns:
        Nothing.

    """

    pathlib.Path().joinpath(
        "results",
        test_configuration["folder"]
    ).mkdir(parents=False, exist_ok=True)

    run_path = pathlib.Path().joinpath(
        "runs",
        test_configuration["folder"]
    )
    run_path.mkdir(parents=False, exist_ok=True)

    pathlib.Path("in").mkdir(parents=False, exist_ok=True)
    pathlib.Path("out").mkdir(parents=False, exist_ok=True)
    pathlib.Path("outb").mkdir(parents=False, exist_ok=True)

def copy_input_test_files(
        test_configuration: dict,
        executable_path: pathlib.Path
) -> None:
    """
    Copies all test input files.

    Copies the pre-compiled executable from the provided executable path to
    the current folder and all the executable input files required by the
    test configuration that were previously fetched.

    Args:
        test_configuration: Dictionary with all test configuration.
        executable_path: Path to the precompiled executable.

    Returns:
        Nothing.

    """

    if pathlib.Path(executable_path.name).exists():
        log.info("Deleting previous executable...")
        pathlib.Path(executable_path.name).unlink()

    log.info("Copying new executable...")
    shutil.copy(
        executable_path,
        executable_path.name
    )

    for input_file in test_configuration["input_files"]:

        shutil.copy(
            pathlib.Path().joinpath(
                "gold_runs",
                test_configuration["folder"],
                input_file
            ),
            pathlib.Path().joinpath(
                "in",
                input_file
            )
        )

def copy_output_test_files(
        test_configuration: dict
):
    """
    Copies output test files.

    Copies the important test output files (temperature map and magnetic field
    map) from the temporary run folder to the test destination run folder.

    Args:
        test_configuration: Dictionary with all test configuration.

    Returns:
        Nothing.

    """

    shutil.copy2(
        pathlib.Path().joinpath(
            "out",
            test_configuration["temperature_map"]
        ),
        pathlib.Path().joinpath(
            "runs",
            test_configuration["folder"],
            test_configuration["temperature_map"]
        )
    )

    shutil.copy2(
        pathlib.Path().joinpath(
            "outb",
            test_configuration["magnetic_field"]
        ),
        pathlib.Path().joinpath(
            "runs",
            test_configuration["folder"],
            test_configuration["magnetic_field"]
        )
    )

    if "magnetic_radial_profiles" in test_configuration:

        for profile in test_configuration["magnetic_radial_profiles"]:

            shutil.copy2(
                pathlib.Path().joinpath(
                    "outb",
                    profile
                ),
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    if "magnetic_meridional_profiles" in test_configuration:

        for profile in test_configuration["magnetic_meridional_profiles"]:

            shutil.copy2(
                pathlib.Path().joinpath(
                    "outb",
                    profile
                ),
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    if "temperature_radial_profiles" in test_configuration:

        for profile in test_configuration["temperature_radial_profiles"]:

            shutil.copy2(
                pathlib.Path().joinpath(
                    "out",
                    profile
                ),
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    if "temperature_meridional_profiles" in test_configuration:

        for profile in test_configuration["temperature_meridional_profiles"]:

            shutil.copy2(
                pathlib.Path().joinpath(
                    "out",
                    profile
                ),
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

def run_test(
        test_name: str,
        test_configuration: dict,
        executable_path: pathlib.Path
) -> None:
    """
    Test running routine.

    This is the main routine for running a particular test. It checks out all
    the required files for the test as specified in its configuration,
    computes all the statistics, and generates comparison plots.

    Args:
        test_name: Name for the test.
        test_configuration: Dictionary with the test configuration.

    Returns:
        Nothing.

    """

    log.info("\n")
    log.info("*** Running test {} *************************".format(test_name))
    log.info("Description: {}\n".format(test_configuration["comments"]))

    # Preconditions, check that the needed files are present.
    checkout_test_files(test_configuration)

    # Create folders for execution.
    create_test_execution_folders(test_configuration)

    # Copy executable and input files.
    copy_input_test_files(test_configuration, executable_path)

    # Launch executable.
    log.info("Running {}!".format(executable_path.name))

    process = subprocess.Popen(
        "./" + executable_path.name,
        cwd="."
    )
    process.wait()

    log.info("Process finished...")

    # Copy output files to appropriate locations.
    copy_output_test_files(test_configuration)

    # Read data from files ----------------------------------------------------
    test_snapshot = test_configuration["test_snapshot"]

    (
        gold_bfield_polar_coordinates,
        gold_bfield_bphi,
        _,
        _
    ) = read_magnetic_field(
        pathlib.Path().joinpath(
            "gold_runs",
            test_configuration["folder"],
            test_configuration["magnetic_field"]
        )
    )

    (
        _,
        run_bfield_bphi,
        _,
        _
    ) = read_magnetic_field(
        pathlib.Path().joinpath(
            "runs",
            test_configuration["folder"],
            test_configuration["magnetic_field"]
        )
    )

    (
        gold_tmap_polar_coordinates,
        gold_tmap_temperatures,
        _
    ) = read_temperature(
        pathlib.Path().joinpath(
            "gold_runs",
            test_configuration["folder"],
            test_configuration["temperature_map"]
        )
    )

    (
        _,
        run_tmap_temperatures,
        _
    ) = read_temperature(
        pathlib.Path().joinpath(
            "runs",
            test_configuration["folder"],
            test_configuration["temperature_map"]
        )
    )

    gold_b_radial_profiles = {}
    run_b_radial_profiles = {}

    if "magnetic_radial_profiles" in test_configuration:

        for profile in test_configuration["magnetic_radial_profiles"]:

            gold_b_radial_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                )
            )

            run_b_radial_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    gold_b_meridional_profiles = {}
    run_b_meridional_profiles = {}

    if "magnetic_meridional_profiles" in test_configuration:

        for profile in test_configuration["magnetic_meridional_profiles"]:

            gold_b_meridional_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                )
            )

            run_b_meridional_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    gold_t_radial_profiles = {}
    run_t_radial_profiles = {}

    if "temperature_radial_profiles" in test_configuration:

        for profile in test_configuration["temperature_radial_profiles"]:

            gold_t_radial_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                )
            )

            run_t_radial_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    gold_t_meridional_profiles = {}
    run_t_meridional_profiles = {}

    if "temperature_meridional_profiles" in test_configuration:

        for profile in test_configuration["temperature_meridional_profiles"]:

            gold_t_meridional_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "gold_runs",
                    test_configuration["folder"],
                    profile
                )
            )

            run_t_meridional_profiles[profile], _, _ = read_profile(
                pathlib.Path().joinpath(
                    "runs",
                    test_configuration["folder"],
                    profile
                )
            )

    # Read grid size from input file, the first input files should contain
    # the input.dat file.
    _, lmax = read_grid_size(
        pathlib.Path().joinpath(
            "in",
            test_configuration["input_files"][0]
        )
    )

    # Polar coordinates are verbose in the snapshot in the sense that there is
    # one line to specify the temperature value for each pair of coordinates.
    # This means that polar coordinates are repeated a lot in each snapshot and
    # we actually only need the ranges and its steps but not each pair of polar
    # coordinates.
    angles_tmap = gold_tmap_polar_coordinates[0, ::(lmax - 1), 0]
    radial_tmap = gold_tmap_polar_coordinates[0, 0:(lmax - 1), 1]
    # The same happens for the magnetic field but this time the grid is
    # complete and does not lack the interface components.
    angles_bfield = gold_bfield_polar_coordinates[0, ::(lmax * 2), 0]
    radial_bfield = gold_bfield_polar_coordinates[0, 0:(lmax * 2), 1]

    # Select snapshots to compare ---------------------------------------------

    gold_bfield = gold_bfield_bphi[test_snapshot]
    gold_tmap = gold_tmap_temperatures[test_snapshot]

    run_bfield = run_bfield_bphi[test_snapshot]
    run_tmap = run_tmap_temperatures[test_snapshot]

    # Preprocess data ---------------------------------------------------------

    temperature_factor = 1e8
    gold_tmap /= temperature_factor
    run_tmap /= temperature_factor

    # Compute stats -----------------------------------------------------------

    error_epsilon = 1e-8

    bfield_abs_error = np.abs(gold_bfield - run_bfield)
    bfield_relative_error = np.divide(
        bfield_abs_error,
        np.maximum(np.abs(gold_bfield), error_epsilon),
    )
    bfield_mean_abs_error = np.mean(bfield_abs_error)

    bfield_relative_l2 = relative_l2_norm(gold_bfield, run_bfield)

    tmap_abs_error = np.abs(gold_tmap - run_tmap)
    tmap_relative_error = np.divide(
        tmap_abs_error,
        np.maximum(np.abs(gold_tmap), error_epsilon),
    )
    tmap_mean_abs_error = np.mean(tmap_abs_error)

    tmap_relative_l2 = relative_l2_norm(gold_tmap, run_tmap)

    # Generate report statistics ----------------------------------------------

    log.info(
        "Magnetic field relative L2 norm: {}".format(
            bfield_relative_l2
        )
    )

    log.info(
        "Magnetic field mean absolute error: {} [G] x 10^12".format(
            bfield_mean_abs_error
        )
    )

    log.info(
        "Temperature map relative L2 norm {}".format(
            tmap_relative_l2
        )
    )

    log.info(
        "Temperature map mean absolute error: {} [K] x 10^8".format(
            tmap_mean_abs_error
        )
    )

    # Plot and save profile comparisons ---------------------------------------

    if "magnetic_radial_profiles" in test_configuration:

        plot_profile_errors(
            gold_b_radial_profiles,
            run_b_radial_profiles,
            None,
            None,
            pathlib.Path().joinpath(
                "results",
                test_configuration["folder"],
                "magnetic_radial_profiles.png"
            )
        )

        for profile in test_configuration["magnetic_radial_profiles"]:

            plot_profile_error(
                gold_b_radial_profiles[profile],
                run_b_radial_profiles[profile],
                pathlib.Path().joinpath(
                    "results",
                    test_configuration["folder"],
                    profile + ".png"
                )
            )

    if "magnetic_meridional_profiles" in test_configuration:

        plot_profile_errors(
            gold_b_meridional_profiles,
            run_b_meridional_profiles,
            None,
            None,
            pathlib.Path().joinpath(
                "results",
                test_configuration["folder"],
                "magnetic_meridional_profiles.png"
            )
        )

        for profile in test_configuration["magnetic_meridional_profiles"]:

            plot_profile_error(
                gold_b_meridional_profiles[profile],
                run_b_meridional_profiles[profile],
                pathlib.Path().joinpath(
                    "results",
                    test_configuration["folder"],
                    profile + ".png"
                )
            )

    if "temperature_radial_profiles" in test_configuration:

        plot_profile_errors(
            gold_t_radial_profiles,
            run_t_radial_profiles,
            None,
            None,
            pathlib.Path().joinpath(
                "results",
                test_configuration["folder"],
                "temperature_radial_profiles.png"
            )
        )

        for profile in test_configuration["temperature_radial_profiles"]:

            plot_profile_error(
                gold_t_radial_profiles[profile],
                run_t_radial_profiles[profile],
                pathlib.Path().joinpath(
                    "results",
                    test_configuration["folder"],
                    profile + ".png"
                )
            )

    if "temperature_meridional_profiles" in test_configuration:

        plot_profile_errors(
            gold_t_meridional_profiles,
            run_t_meridional_profiles,
            None,
            None,
            pathlib.Path().joinpath(
                "results",
                test_configuration["folder"],
                "temperature_meridional_profiles.png"
            )
        )

        for profile in test_configuration["temperature_meridional_profiles"]:

            plot_profile_error(
                gold_t_meridional_profiles[profile],
                run_t_meridional_profiles[profile],
                pathlib.Path().joinpath(
                    "results",
                    test_configuration["folder"],
                    profile + ".png"
                )
            )

    # Plot and save eror maps -------------------------------------------------

    plot_error(
        angles_tmap,
        radial_tmap,
        tmap_abs_error,
        tmap_relative_error,
        "T ~ [K] x 10^8",
        pathlib.Path().joinpath(
            "results",
            test_configuration["folder"],
            "temperature_error.png"
        )
    )

    log.info(
        "Temperature plot generated in results/{}".format(
            test_configuration["folder"]
        )
    )

    plot_error(
        angles_bfield,
        radial_bfield,
        bfield_abs_error,
        bfield_relative_error,
        "B ~ [G] x 10^{12}",
        pathlib.Path().joinpath(
            "results",
            test_configuration["folder"],
            "magnetic_error.png"
        )
    )

    log.info(
        "Magnetic field plot generated in results/{}".format(
            test_configuration["folder"]
        )
    )


def run_tests(
        regression_config: dict,
        executable_path: pathlib.Path
) -> None:
    """
    Run the whole test suite as specified in the configuration.

    Loops over the test configuration dictionary extracted from the JSON config
    file and runs each entry one by one.

    Args:
        regression_config: JSON dictionary with the test configuration

    Returns:
        Nothing.

    """

    log.info("Running regression tests...")

    for test_name, test_configuration in regression_config.items():
        run_test(test_name, test_configuration, executable_path)


def cleanup(
        executable_path: pathlib.Path,
) -> None:
    """
    Clean up all files.

    Args:
        executable_path: Path to the executable.

    Returns:
        Nothing.
    """

    # Remove executable.
    if pathlib.Path(executable_path.name).exists():
        pathlib.Path(executable_path.name).unlink()

    # Remove data from folders.
    for root, dirs, _ in os.walk("results"):
        for directory in dirs:
            shutil.rmtree(os.path.join(root, directory))

    for root, dirs, _ in os.walk("gold_runs"):
        for directory in dirs:
            shutil.rmtree(os.path.join(root, directory))

    for root, dirs, _ in os.walk("runs"):
        for directory in dirs:
            shutil.rmtree(os.path.join(root, directory))

    # Remove temporary data.
    if pathlib.Path("in").exists():
        shutil.rmtree(pathlib.Path("in"))
    if pathlib.Path("out").exists():
        shutil.rmtree(pathlib.Path("out"))
    if pathlib.Path("outb").exists():
        shutil.rmtree(pathlib.Path("outb"))


if __name__ == "__main__":

    # Parse arguments ---------------------------------------------------------

    parser = argparse.ArgumentParser(description="Parameters")
    parser.add_argument(
        "--log_file_path",
        nargs="?",
        type=str,
        default="results/log.txt",
        help="Log file path"
    )
    parser.add_argument(
        "--regression_configuration_filename",
        nargs="?",
        type=str,
        default="regression_configuration.json",
        help="Regression configuration filename"
    )
    parser.add_argument(
        "--cleanup",
        nargs="?",
        type=bool,
        default=False,
        help="Cleanup all files before running"
    )
    args = parser.parse_args()

    LOG_FORMAT = "%(asctime)s [%(levelname)-5.5s]  %(message)s"
    LOG_FILE_PATH = pathlib.Path(args.log_file_path)
    REGRESSION_CONFIGURATION_FILENAME = args.regression_configuration_filename
    DO_CLEANUP = args.cleanup
    EXECUTABLE_PATH = pathlib.Path("../../build/bin/mt2d")

    # Setup logger ------------------------------------------------------------

    # Delete existing log file if it exists.
    if LOG_FILE_PATH.exists():
        LOG_FILE_PATH.unlink()

    # Configure standard logger to console.
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    # Configure logging to file.
    log_formatter = logging.Formatter(LOG_FORMAT)
    log_file_handler = logging.FileHandler(LOG_FILE_PATH)
    log_file_handler.setFormatter(log_formatter)
    log.addHandler(log_file_handler)

    # Configure plotting ------------------------------------------------------

    configure_plotting()

    # Data reading ------------------------------------------------------------

    REGRESSION_CONFIGURATION = read_regression_configuration(
        REGRESSION_CONFIGURATION_FILENAME
    )

    if DO_CLEANUP:
        cleanup(EXECUTABLE_PATH)

    # Test suite execution ----------------------------------------------------

    run_tests(REGRESSION_CONFIGURATION, EXECUTABLE_PATH)
