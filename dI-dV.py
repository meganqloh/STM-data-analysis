# Import necessary libraries
import sys
from PyQt5 import QtWidgets
from pathvalidate import sanitize_filename
from time import ctime
from numpy import savetxt
from math import dist
from PyQt5.QtCore import Qt

# Import user-defined library
from omiclib import Data3D, dIdVGraph, FFTGraph, FFTLineCutGraph, SliderAndFloatInput, RangeSliderAndFloatInputs, PointSetter, format_file_names, angle

# Define PromptDialog class for selecting options and files
class PromptDialog(QtWidgets.QDialog):
    def __init__(self, options, parent=None):
        super().__init__(parent)
        
        self.file_names = []
        self.selection = None
        self.setWindowTitle('dI/dV Maps Init')
        
        # Create QPushButton for file selection
        file_label = QtWidgets.QLabel('dI/dV file(s) to analyze (.Aux2(V)_mtrx or .asc):')
        file_button = QtWidgets.QPushButton('Select file(s)')
        file_button.clicked.connect(self.selectFile)
        
        # Create QComboBox for the options
        question_label = QtWidgets.QLabel('Function to use:')
        self.option_combobox = QtWidgets.QComboBox()
        self.option_combobox.addItems(options)
        self.selected = QtWidgets.QLabel()
        
        # Create QPushButton for confirmation
        confirm_button = QtWidgets.QPushButton('Confirm')
        confirm_button.clicked.connect(self.confirmChoice)
        
        self.status = QtWidgets.QStatusBar()
        
        # Create QVBoxLayout to arrange the widgets vertically
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(file_label)
        layout.addWidget(file_button)
        layout.addWidget(self.selected)
        layout.addWidget(question_label)
        layout.addWidget(self.option_combobox)
        layout.addSpacing(50)
        layout.addWidget(confirm_button)
        layout.addWidget(self.status)
        
        self.setLayout(layout)
        
    def confirmChoice(self):
        # Get the selected option from the QComboBox
        self.selection = self.option_combobox.currentText()
        try:
            # Check if files are selected and their types are correct
            if len(self.file_names) > 0:
                if not (len(self.file_names) == 1 and self.file_names[0].endswith('.Aux2(V)_mtrx')):
                    self.file_names = [x for x in self.file_names if x.endswith('.asc')]
                    if len(self.file_names) > 0:
                        self.close()
                else:
                    self.close()
            else:
                self.status.showMessage("Incorrect file type / no file selected.", int(1e4))
        except Exception as e:
            self.status.showMessage(str(e), int(1e4))
            
    def selectFile(self):
        # Open a file dialog to select files
        self.file_names, _ = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select .Aux2(V)_mtrx/.asc file(s)')
        self.selected.setText(format_file_names(self.file_names))

# Define MainWindow class for the main application window
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, prompt_selection, data):
        super().__init__()
        
        self.setWindowTitle('dI/dV Map Interface')
        selection = [r'Only dI/dV maps', r'dI/dV maps and FFT', r'dI/dV FFT line cuts']
        self.op = selection.index(prompt_selection)
        self.V = data.V
        self.energy_dict = data.energy_dict
        init_energy = min(data.V[0], data.V[-1])
            
        # Initialize all widgets for each function
        if len(self.V) > 1:
            self.energy_slider = SliderAndFloatInput('Choose energy', data.V, self)
        self.didv = dIdVGraph(self.energy_dict.get(init_energy), init_energy)
        if self.op == 1:
            self.fft = FFTGraph(self.energy_dict.get(init_energy))
            self.fft_slider = RangeSliderAndFloatInputs('Adjust colormap', 0, self.fft.getMax(), self)
        if self.op == 2:
            self.fftlc = FFTLineCutGraph(data, init_energy, self)
            self.fft_slider = RangeSliderAndFloatInputs('Adjust colormap of line cut/FFT plots', 0, self.fftlc.getMax(), self)
            self.qpi_slider = RangeSliderAndFloatInputs('Adjust colormap of QPI plot', 0, self.fftlc.getMax(), self, 2)
            # Formatting and initializing line cut point setter and text with angle and distance
            self.lc_setpt = PointSetter(self.fftlc)
            self.lc_coords = QtWidgets.QLabel()
            pt_layout = QtWidgets.QVBoxLayout()
            pt_layout.addWidget(self.lc_setpt)
            pt_layout.addWidget(self.lc_coords)
            pt = QtWidgets.QWidget()
            pt.setLayout(pt_layout)
            self.updateCoords()
        
        # Format main widgets
        main_layout = QtWidgets.QHBoxLayout()
        main_layout.addWidget(self.didv)
        if self.op == 1:
            main_layout.addWidget(self.fft)
        if self.op == 2:
            main_layout.addWidget(self.fftlc)
            main_layout.setStretchFactor(self.didv, 1)
            main_layout.setStretchFactor(self.fftlc, 2)
        main = QtWidgets.QWidget()
        main.setLayout(main_layout)
        self.setCentralWidget(main)
        
        # Format dock widgets
        dock = QtWidgets.QDockWidget('Make selections', self)
        dock.setFloating(False)
        dock_splitter = QtWidgets.QSplitter(self)
        if len(self.V) > 1:
            dock_splitter.addWidget(self.energy_slider)
        if self.op >= 1:
            dock_splitter.addWidget(self.fft_slider)
            if self.op == 2:
                dock_splitter.addWidget(pt)
                dock_splitter.addWidget(self.qpi_slider)
        dock.setWidget(dock_splitter)
        self.addDockWidget(Qt.TopDockWidgetArea, dock)
        
        self.file_identifier = QtWidgets.QLineEdit()
        self.file_identifier.setFixedWidth(190)
        self.save_didv = QtWidgets.QPushButton('Save all dI/dV maps', self)
        self.save_didv.clicked.connect(lambda: self.saveData('all'))
        self.save_didv.setEnabled(True)
        if self.op == 2:
            self.save_lc_qpi = QtWidgets.QPushButton('Save line cut / QPI', self)
            self.save_lc_qpi.clicked.connect(lambda: self.saveData('LC QPI'))
            self.save_lc_qpi.setEnabled(True)
        
        fi_toolbar = self.addToolBar('Choose file identifier')
        fi_toolbar.addWidget(QtWidgets.QLabel('File identifier: /selected/path/'))
        fi_toolbar.addWidget(self.file_identifier)
        label = '_dIdVmap_#.##V.asc'
        if self.op == 2:
            label += ' OR _line_cut_#.##V.asc OR _qpi.asc'
        label += '   '
        fi_toolbar.addWidget(QtWidgets.QLabel(label))
        self.addToolBarBreak()
        save_toolbar = self.addToolBar('Save data')
        save_toolbar.addWidget(self.save_didv)
        if self.op == 2:
            save_toolbar.addWidget(self.save_lc_qpi)
            
        self.status_bar = self.statusBar()
    
    def saveData(self, button):
        # Save the data based on the selected option
        pre_text = sanitize_filename(self.file_identifier.text())
        dialog = QtWidgets.QFileDialog()
        folder = dialog.getExistingDirectory(None, "Select folder to save inside")
        if folder is not None:
            save_message = 'None'
            if button == 'all':
                for E, data in self.energy_dict.items():
                    header = f'File format = ASCII\nx-pixels = {data.shape[1]}\ny-pixels = {data.shape[0]}\nEnergy = {E}\nEnergy units = V\ndI/dV maps\n\nStart of data:'
                    if len(pre_text) > 0:
                        filename = f'{pre_text}_dIdVmap_{E}V.asc'
                    else:
                        filename = f'dIdVmap_{E}V.asc'
                    savetxt(f'{folder}/{filename}', data, delimiter='\t', header=header, fmt='%.6f')
                save_message = 'All dI/dV maps'
            if button == 'LC QPI':
                self.fftlc.saveData(folder, pre_text)
                save_message = 'Line cut + QPI'
            save_message += f' saved at {ctime()}'
            self.status_bar.showMessage(save_message, int(1e4))
                
    
    def updatePlots(self, value):
        # Update the plots based on the selected energy value
        self.didv.updateFigure(self.energy_dict.get(value), value)
        if self.op == 1:
            self.fft.updateFigure(self.energy_dict.get(value))
        if self.op == 2:
            self.fftlc.updateFigure(value)
            
    def updateCoords(self):
        # Update the QLabel showing angle and distance of line cut
        pt1 = self.fftlc.list_points[0].point.center
        pt2 = self.fftlc.list_points[1].point.center
        self.lc_coords.setText(f'Angle: {round(angle(pt1[0], pt1[1], pt2[0], pt2[1]), 1)}\u00b0, Distance: {round(dist(pt1, pt2), 2)}')
            
    def updateCmap(self, lower, upper, slider_id):
        # Update the colormap based on slider values
        if slider_id == 1:
            if self.op == 1:
                self.fft.updateClim(lower, upper)
            if self.op == 2:
                self.fftlc.updateClimFFTLC(lower, upper)
                self.fftlc.draw_idle()
        if slider_id == 2:
            if self.op == 2:
                self.fftlc.updateClimQPI(lower, upper)
                self.fftlc.draw_idle()
    
    def getClim(self, slider_id):
        # Get the current colormap limits for the selected slider
        if slider_id == 1:
            return self.fft_slider.getRange()
        if slider_id == 2:
            return self.qpi_slider.getRange()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    
    # Prompt for the option
    options = [r'Only dI/dV maps', r'dI/dV maps and FFT', r'dI/dV FFT line cuts']
    init_prompt = PromptDialog(options)
    init_prompt.exec_()
    data = Data3D(init_prompt.file_names)
    if not data.valid:
        print('Data not compatible.')
        sys.exit()
    window = MainWindow(init_prompt.selection, data)
    window.show()
    
    sys.exit(app.exec_())
