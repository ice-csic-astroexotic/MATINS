#! /usr/bin/python3

"""
MagPyGraph.

This program replicates the functionality of `ygraph` for the MAGNESIA project.

    Authors:
        Alberto Garcia Garcia (garciagarcia@ice.csic.es)

MIT License

Copyright (c) MAGNESIA (ICE-CSIC) 2020

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse
import logging
import pathlib
import sys
import typing

import numpy as np
from matplotlib.pyplot import cm
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtWidgets import (
    QWidget,
    QHBoxLayout,
    QPushButton,
    QGridLayout,
    QLabel,
    QLineEdit,
    QCheckBox
)
import pyqtgraph as pg
# from pyqtgraph.Qt import QtGui

MAX_PROFILES = 8
BACKGROUND_COLOR = (255, 255, 255) 
LINE_COLOR = (255, 255, 255)
LINE_WIDTH = 2.0
LINE_STYLE = QtCore.Qt.SolidLine
COLOR_MAP_ITERATOR = iter(cm.hsv(np.linspace(0, 1, MAX_PROFILES)))
SYMBOL = 'x'
SYMBOL_SIZE = 4
BUTTON_WIDTH = 80 #128
BUTTON_HEIGHT = 64


logging.basicConfig(level=logging.DEBUG)


def read_profile(
        filename: str,
) -> typing.Tuple[np.array, np.array, np.array]:
    """
    Auxiliary function for reading a ygraph profile `.yg`.

        Args:
            filename: path to the profile to load.

        Returns:
            A triplet of numpy arrays with the values and their coordinates for
            each snapshot in the profile together with the timestamps for them.
    """

    coordinates: list = []
    values: list = []
    times: list = []

    with open(filename, "r") as _f:

        current_coordinates: list = []
        current_values: list = []

        # Skip the first line which contains the label.
        next(_f)

        for line in _f:

            split = line.strip().split()

            # If the line is a blank line, we store.
            if len(split) == 0:

                coordinates.append(current_coordinates)
                values.append(current_values)

                current_coordinates = []
                current_values = []

            # If the line is a "Label" we skip it.
            elif "\"Label" in split[0]:

                continue

            # If the line is a time tag, store it.
            elif "\"Time=" in split[0]:

                times.append(float(split[1]))

            # Otherwise, the line is a coordinate-value pair.
            else:

                current_coordinates.append(float(split[0]))
                current_values.append(float(split[1]))

    coordinates = np.array(coordinates)
    values = np.array(values)
    times = np.array(times)

    return times, coordinates, values


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # UI ------------------------------------------------------------------

        # Graph.
        self.graphWidget = pg.PlotWidget()
        self.graphWidget.setBackground(BACKGROUND_COLOR)
        self.graphWidget.showGrid(x=True, y=True)
        self.graphWidget.setLogMode(False, False)
        self.legend = self.graphWidget.addLegend()
        # Buttons.
        self.button_next = QPushButton("Next")
        self.button_next.setFixedWidth(BUTTON_WIDTH)
        self.button_next.setFixedHeight(BUTTON_HEIGHT)
        self.button_next.clicked.connect(self.next_plot_data)
        self.button_end = QPushButton("End")
        self.button_end.setFixedWidth(BUTTON_WIDTH)
        self.button_end.setFixedHeight(BUTTON_HEIGHT)
        self.button_end.clicked.connect(self.go_to_end)
        self.button_prev = QPushButton("Previous")
        self.button_prev.setFixedWidth(BUTTON_WIDTH)
        self.button_prev.setFixedHeight(BUTTON_HEIGHT)
        self.button_prev.clicked.connect(self.prev_plot_data)
        self.button_beg = QPushButton("Start")
        self.button_beg.setFixedWidth(BUTTON_WIDTH)
        self.button_beg.setFixedHeight(BUTTON_HEIGHT)
        self.button_beg.clicked.connect(self.go_to_beginning)
        self.button_play = QPushButton("Play")
        self.button_play.setFixedWidth(BUTTON_WIDTH)
        self.button_play.setFixedHeight(BUTTON_HEIGHT)
        self.button_play.clicked.connect(self.start_timer)
        self.button_stop = QPushButton("Stop")
        self.button_stop.setFixedWidth(BUTTON_WIDTH)
        self.button_stop.setFixedHeight(BUTTON_HEIGHT)
        self.button_stop.clicked.connect(self.stop_timer)
        self.button_rescale = QPushButton("Rescale Y to Timestep")
        self.button_rescale.setFixedWidth(BUTTON_WIDTH)
        self.button_rescale.setFixedHeight(BUTTON_HEIGHT)
        self.button_rescale.clicked.connect(self.rescale_y_to_timestep)
        self.button_rescale_seq = QPushButton("Rescale Y to Sequence")
        self.button_rescale_seq.setFixedWidth(BUTTON_WIDTH)
        self.button_rescale_seq.setFixedHeight(BUTTON_HEIGHT)
        self.button_rescale_seq.clicked.connect(self.rescale_y_to_sequence)
        self.button_rescale_x = QPushButton("Rescale X to All")
        self.button_rescale_x.setFixedWidth(BUTTON_WIDTH)
        self.button_rescale_x.setFixedHeight(BUTTON_HEIGHT)
        self.button_rescale_x.clicked.connect(self.rescale_x_to_coord_range)
        self.button_recenter = QPushButton("Recenter")
        self.button_recenter.setFixedWidth(BUTTON_WIDTH)
        self.button_recenter.setFixedHeight(BUTTON_HEIGHT)
        self.button_recenter.clicked.connect(self.recenter)
        self.button_goto = QPushButton("Go to Frame")
        self.button_goto.setFixedWidth(BUTTON_WIDTH)
        self.button_goto.setFixedHeight(BUTTON_HEIGHT)
        self.button_goto.clicked.connect(self.go_to_specified_frame)
        # Labels.
        self.label_timestep = QLabel()
        self.label_timestep.setText("Timestep: 0")
        self.label_time = QLabel()
        self.label_time.setText("Time: 0")
        # Textbox.
        self.text_frame_selector = QLineEdit()
        self.text_frame_selector.setFixedWidth(BUTTON_WIDTH)
        self.text_frame_selector.setFixedHeight(BUTTON_HEIGHT)
        # Checkbox.
        self.checkbox_logscaley = QCheckBox("Log Scale Y-axis")
        self.checkbox_logscaley.stateChanged.connect(
            lambda: self.graphWidget.setLogMode(
                y=self.checkbox_logscaley.isChecked()
            )
        )
        self.checkbox_logscalex = QCheckBox("Log Scale X-axis")
        self.checkbox_logscalex.stateChanged.connect(
            lambda: self.graphWidget.setLogMode(
                x=self.checkbox_logscalex.isChecked()
            )
        )

        # Grid.
        self.grid = QGridLayout()

        self.button_box = QHBoxLayout()
        self.button_box.addWidget(self.button_goto)
        self.button_box.addWidget(self.text_frame_selector)
        self.button_box.addWidget(self.button_beg)
        self.button_box.addWidget(self.button_prev)
        self.button_box.addWidget(self.button_next)
        self.button_box.addWidget(self.button_end)
        self.button_box.addWidget(self.button_stop)
        self.button_box.addWidget(self.button_play)
        self.button_box.addWidget(self.label_timestep)
        self.button_box.addWidget(self.label_time)
        self.button_box.addWidget(self.checkbox_logscaley)
        self.button_box.addWidget(self.checkbox_logscalex)
        self.button_box.addWidget(self.button_rescale)
        self.button_box.addWidget(self.button_rescale_seq)
        self.button_box.addWidget(self.button_rescale_x)
        self.button_box.addWidget(self.button_recenter)

        self.grid.addLayout(self.button_box, 0, 0)
        self.grid.addWidget(self.graphWidget, 1, 0)

        # Widget.
        self.widget = QWidget()
        self.widget.setLayout(self.grid)
        self.setCentralWidget(self.widget)

        # Change legend font size.
        legendLabelStyle = {
            'color':'#FFF',
            'size': '24pt',
            'bold': True,
            'italic': False
        }

        # Set legend and ticks font size.
        self.legend.setScale(1.4)
        #font = QtGui.QFont()
        #font.setPixelSize(24)
        #self.graphWidget.getAxis("bottom").setTickFont(font)
        #self.graphWidget.getAxis("bottom").setStyle(tickTextOffset = 20)
        #self.graphWidget.getAxis("left").setTickFont(font)
        #self.graphWidget.getAxis("left").setStyle(tickTextOffset = 20)

        # Plotting data -------------------------------------------------------
        self.plots: dict = {}
        self.profiles: dict = {}
        self.max_timestep: int = sys.maxsize
        self.timestep: int = 0
        self.times: list = []
        self.coordinates: list = []

        # Timer ---------------------------------------------------------------
        self.timer = QtCore.QTimer()
        self.timer.setInterval(50)
        self.timer.timeout.connect(self.next_plot_data)

    def go_to_end(self) -> None:
        """ Stops the timer and sets the current frame to the last one. """
        self.stop_timer()
        self.go_to_frame(self.max_timestep)

    def go_to_beginning(self) -> None:
        """ Stops the timer and sets the current frame to the first one. """
        self.stop_timer()
        self.go_to_frame(0)

    def go_to_frame(self, frame: int) -> None:
        """ Sets the current frame to the specified one for all plots. """

        self.timestep = frame

        for profile_name, data in self.profiles.items():
            x = data["x"][self.timestep]
            y = data["y"][self.timestep]
            self.plots[profile_name].setData(x, y)

        self.label_timestep.setText(f"Timestep: {self.timestep}")
        self.label_time.setText(f"Time: {self.times[self.timestep]:.0f}")

    def load_profile(
        self,
        profile_path: pathlib.Path
    ) -> typing.Tuple[str, int]:

        times, x_values, y_values = read_profile(profile_path)
        profile_name = profile_path.name

        self.profiles[profile_name] = {
            "time": times,
            "x": x_values,
            "y": y_values
        }

        profile_min = np.min(y_values)
        profile_max = np.max(y_values)

        logging.info(f"Loaded profile {profile_name}")
        logging.info(f"Profile range is [{profile_min}, {profile_max}]...")

        return profile_name, len(times) - 1

    def load_profiles(self, profile_list: typing.List[str]) -> None:

        for profile in profile_list:

            profile_path = pathlib.Path(profile)
            # TODO: verify path exists and warn accordingly.
            profile_name, timesteps = self.load_profile(profile_path)
            self.max_timestep = min(self.max_timestep, timesteps)

            self.plots[profile_name] = self.plot_line(
                profile_name,
                self.profiles[profile_name]["x"][self.timestep],
                self.profiles[profile_name]["y"][self.timestep],
                next(COLOR_MAP_ITERATOR) * 255.0,
                SYMBOL,
                LINE_WIDTH,
                LINE_STYLE
            )

        logging.info(f"Maximum common timestep is {self.max_timestep}")

    def verify_times(self):
        """ Verifies that all loaded profiles contain the same timestamps. """
        times: list = []

        for profile_name, data in self.profiles.items():
            times.append(data["time"][:self.max_timestep+1])

        if not all(np.array_equal(x, times[0]) for x in times):
            logging.error("Profiles are not synchronized...")
            sys.exit()

        self.times = times[0]

    def verify_coordinates(self):
        """ Verifies that all loaded profiles contain the same coordinates. """

        coordinates: list = []

        for profile_name, data in self.profiles.items():
            coordinates.append(data["x"][0])

        if not all(np.array_equal(x, coordinates[0]) for x in coordinates):
            logging.error("Profiles coordinates do not match...")
            sys.exit()

        self.coordinates = coordinates[0]

    def plot_line(self, name, x, y, color, symbol, width, style):
        pen = pg.mkPen(color=(color), width=width, style=style)
        brush = pg.mkBrush(color=(color))
        return self.graphWidget.plot(
            x,
            y,
            name=name,
            pen=pen,
            symbol=symbol,
            symbolPen=pen,
            symbolBrush=brush,
            symbolSize=SYMBOL_SIZE
        )

    # Callbacks ---------------------------------------------------------------

    def go_to_specified_frame(self) -> None:
        """ Callback function for the frame GOTO button, reads the specified
            frame from the textbox and calls the routine to set it.
        """
        self.stop_timer()
        frame = int(self.text_frame_selector.text())
        self.go_to_frame(frame)

    def recenter(self) -> None:
        """ Callback to rescale both axes. """
        self.rescale_x_to_coord_range()
        self.rescale_y_to_sequence()

    def rescale_x_to_coord_range(self) -> None:
        """ Callback to rescale x-axis to the min/max coordinate range. """
        min_value: float = min(self.coordinates)
        max_value: float = max(self.coordinates)

        self.graphWidget.setXRange(min_value, max_value, padding=0)

        logging.info(f"Rescaled range to x = [{min_value}, {max_value}]")

    def rescale_y_to_timestep(self) -> None:
        """ Callback to rescale y-axis to the current frame min/max range. """
        min_value: float = float("inf")
        max_value: float = float("-inf")

        for name, data in self.profiles.items():
            min_value = min(np.min(data["y"][self.timestep]), min_value)
            max_value = max(np.max(data["y"][self.timestep]), max_value)

        self.graphWidget.setYRange(min_value, max_value, padding=0)

        logging.info(f"Rescaled range to y = [{min_value}, {max_value}]")

    def rescale_y_to_sequence(self) -> None:
        """ Callback to rescale y-axis to the whole sequence min/max range. """
        min_value: float = float("inf")
        max_value: float = float("-inf")

        for name, data in self.profiles.items():
            min_value = min(np.min(data["y"]), min_value)
            max_value = max(np.max(data["y"]), max_value)

        self.graphWidget.setYRange(min_value, max_value, padding=0)

        logging.info(f"Rescaled range to y = [{min_value}, {max_value}]")

    def start_timer(self):
        """ Callback to start playback. """
        self.timer.start()

    def stop_timer(self):
        """ Callback to stop playback. """
        self.timer.stop()

    def next_plot_data(self):
        """ Callback to advance all plots to the next frame. """

        if self.timestep < self.max_timestep:
            self.timestep += 1
            self.go_to_frame(self.timestep)
        else:
            self.stop_timer()

    def prev_plot_data(self):
        """ Callback to rewind all plots to the previous frame. """

        if self.timestep > 0:
            self.timestep -= 1
            self.go_to_frame(self.timestep)
        else:
            self.stop_timer()


def main(args):
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    # Load the specified profiles and verify they all contain the same times
    # and the same coordinates.
    main.load_profiles(args.profiles)
    main.verify_times()
    main.verify_coordinates()
    # Rescale y-axis to the whole sequence min/max values.
    main.rescale_y_to_sequence()
    # Rescale x-axis to the coordinate range.
    main.rescale_x_to_coord_range()
    # Display the plot.
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Parameters")

    parser.add_argument(
        "--profiles",
        nargs="*",
        type=str,
        required=True,
        help="List of profiles")

    args = parser.parse_args()

    main(args)
