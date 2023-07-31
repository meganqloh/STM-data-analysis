#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5 import QtWidgets
import sys
import os

# Personal modules
from omiclib import Data2D, LineCutGraph, PointSetter

class MyWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        # Create the LineCutGraph widget as the central widget of the main window
        self.graph = LineCutGraph()
        self.setCentralWidget(self.graph)

        # Add a button to select the file to analyze
        select_file_button = QtWidgets.QPushButton('Select file to analyze', self)
        select_file_button.clicked.connect(self.selectFile)

        # Add a button to save the plot
        save_button = QtWidgets.QPushButton('Save plot', self)
        save_button.clicked.connect(self.graph.savePlot)
        save_button.setEnabled(False)

        # Add a toolbar for the file operations
        toolbar = self.addToolBar('File Toolbar')
        toolbar.addWidget(select_file_button)
        toolbar.addWidget(save_button)

        # Store the buttons for later use
        self.select_file_button = select_file_button
        self.save_button = save_button

        # Add the navigation toolbar for the LineCutGraph
        self.nav_toolbar = NavigationToolbar(self.graph, self)
        self.addToolBar(self.nav_toolbar)

        # Add the PointSetter widget to set draggable points
        self.addToolBarBreak()
        self.addToolBar('Coordinates').addWidget(PointSetter(self.graph))
        
        self.status = self.statusBar()

    def selectFile(self):
        # Open a file dialog to select the file to analyze
        working_directory = os.getcwd()
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select File', working_directory)
        if filename:
            try:
                # Load and validate the data from the selected file
                data = Data2D(filename)
                if data.valid:
                    self.graph.data = data.data
                    self.graph.x = data.x
                    self.graph.y = data.y
                    # Update the LineCutGraph plots with the data and initial draggable points
                    self.graph.updatePlots([self.graph.x[0] * 2 / 3, self.graph.y[0] * 2 / 3],
                                            [self.graph.x[-1] * 2 / 3, self.graph.y[-1] * 2 / 3])
                    self.save_button.setEnabled(True)
            except:
                # Show an error message if the data is invalid
                self.status.showMessage('Invalid data', int(1e4))
                
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
