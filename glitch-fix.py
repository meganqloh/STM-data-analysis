# Import necessary libraries
from PyQt5.QtWidgets import QApplication, QDialog, QPushButton, QLabel, QVBoxLayout, QHBoxLayout, QStatusBar, QFileDialog, QMainWindow, QDockWidget, QSplitter, QRadioButton, QWidget, QLineEdit
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from numpy import savetxt, loadtxt, roll, fliplr, array
from pathvalidate import sanitize_filename
from time import ctime
import sys

# Import user-defined library
from omiclib import format_file_names, Data3D, SliderAndFloatInput, GlitchFixGraph

# Define the PromptDialog class for selecting files and optional slice
class PromptDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.file_names = ['']
        self.pt_name = ''
        self.setWindowTitle('Glitch Fix Initialization')
        
        # Create QPushButton for file selection
        file_button = QPushButton('Select file(s) to analyze')
        self.selected_files_label = QLabel()
        file_button.clicked.connect(self.selectFile)
        
        # Create QPushButton for selecting optional slice
        pts_button = QPushButton('(Optional) Select file with slice')
        self.selected_pts_label = QLabel()
        pts_button.clicked.connect(self.selectPt)
        
        # Create QPushButton for confirmation
        confirm_button = QPushButton('Confirm')
        confirm_button.clicked.connect(self.confirmChoice)
        
        # Create status bar for error messages
        self.status = QStatusBar()
        
        # Create QVBoxLayout to arrange the widgets vertically
        layout = QVBoxLayout()
        layout.addWidget(QLabel('dI/dV file(s) to analyze (.Aux2(V)_mtrx or .asc):'))
        layout.addWidget(file_button)
        layout.addWidget(self.selected_files_label)
        layout.addWidget(QLabel('(Optional) file with coordinates of slice:'))
        layout.addWidget(pts_button)
        layout.addWidget(self.selected_pts_label)
        layout.addSpacing(20)
        layout.addWidget(confirm_button)
        layout.addWidget(self.status)
        
        self.setLayout(layout)
        
    def confirmChoice(self):
        # Handle confirmation and checking file types
        if len(self.file_names) == 0:
            self.status.showMessage("No file selected.", int(1e4))
            return
        if len(self.file_names) == 1 and self.file_names[0].endswith('.Aux2(V)_mtrx'):
            self.close()
            return
        self.file_names = [x for x in self.file_names if x.endswith('.asc')] # only keep `.asc` files
        if len(self.file_names) > 0:
            self.close()
            return
        self.status.showMessage("Incorrect file type(s).", int(1e4))
            
    def selectFile(self):
        # Select the files for analysis
        self.file_names, _ = QFileDialog.getOpenFileNames(self, 'Select file(s)')
        self.selected_files_label.setText(format_file_names(self.file_names))
        
    def selectPt(self):
        # Select the optional file with coordinates of the slice
        self.pt_name, _ = QFileDialog.getOpenFileName(self, 'Select file')
        self.selected_pts_label.setText(format_file_names([self.pt_name]))

# Define the MainWindow class for the main application window
class MainWindow(QMainWindow):
    def __init__(self, data, pts):
        super().__init__()
        self.setWindowTitle('Glitch Fix')
        self.energy_dict = data.energy_dict
        self.current_energy = min(data.V[0], data.V[-1])
        
        # Create energy slider if more than one dI/dV map is chosen, add to top dock
        if len(data.V) > 1:
            self.energy_slider = SliderAndFloatInput('Choose energy', data.V, self)
            energy_dock = QDockWidget()
            energy_dock.setWidget(self.energy_slider)
            self.addDockWidget(Qt.TopDockWidgetArea, energy_dock)

        # Create main editable plot and set as main widget
        self.graph = GlitchFixGraph(pts, self)
        self.setCentralWidget(self.graph)
        
        # Create matplotlib navigation toolbar
        nav_toolbar = NavigationToolbar(self.graph, self)
        self.addToolBar(nav_toolbar)

        # Create radio buttons to choose between 'Shift' and 'Edit' functions
        self.edit_radio_button = QRadioButton("Edit")
        self.edit_radio_button.clicked.connect(lambda: self.on_radio_button_click("Edit"))
        self.shift_radio_button = QRadioButton("Shift")
        self.shift_radio_button.clicked.connect(lambda: self.on_radio_button_click("Shift"))
        
        # Initialize to 'Edit'
        self.edit_radio_button.setChecked(True)
        self.mode = 'Edit'
        self.accepts_double_click = True

        # Format radio buttons 
        radio_button_layout = QVBoxLayout()
        radio_button_layout.addWidget(QLabel("Mode:"))
        radio_button_layout.addWidget(self.edit_radio_button)
        radio_button_layout.addWidget(self.shift_radio_button)
        radio_button_layout.addStretch()
        radio_widget = QWidget()
        radio_widget.setLayout(radio_button_layout)

        # Create save button
        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.on_save_button_click)
        self.save_button.setVisible(False)

        # Create file identifier textbox to save files under unique file names
        file_identifier_layout = QHBoxLayout()
        file_identifier_layout.addWidget(QLabel('File identifier: /selected/path/'))
        self.file_identifier = QLineEdit()
        file_identifier_layout.addWidget(self.file_identifier)
        file_identifier_layout.addWidget(QLabel('_dIdVmap_#.##V.asc / _glitch_pts.asc   '))

        # Format save button
        save_layout = QVBoxLayout()
        save_layout.addLayout(file_identifier_layout)
        save_layout.addWidget(self.save_button)
        save_widget = QWidget()
        save_widget.setLayout(save_layout)

        # Create bottom dock for radio buttons and save button
        buttons_dock = QDockWidget()
        dock_splitter = QSplitter()
        dock_splitter.addWidget(radio_widget)
        dock_splitter.addWidget(save_widget)
        buttons_dock.setWidget(dock_splitter)
        self.addDockWidget(Qt.BottomDockWidgetArea, buttons_dock)

        # Create status bar for save messages
        self.status_bar = self.statusBar()
        
        # Display first plot on GlitchFixGraph
        self.updatePlots(self.current_energy)

    def updatePlots(self, E):
        # Update data to display and update plots based on active radio button
        self.current_energy = E
        self.graph.update_data(self.current_energy)
        self.on_radio_button_click(self.mode)

    def on_radio_button_click(self, label):
        self.mode = label
        if label == 'Shift': # Turns off graph interactivity, shifts all data, allows saving
            self.accepts_double_click = False
            self.graph.shift_data()
            self.save_button.setVisible(True)
        if label == 'Edit': # Turns on graph interactivity
            self.accepts_double_click = True
            self.save_button.setVisible(False)
        self.graph.update_plot()
    
    def on_save_button_click(self, event):
        # Save the shifted data to files
        sanitized_identifier = sanitize_filename(self.file_identifier.text()) # Ensures file identifier can be used in filename
        
        # Prompt user to select folder to save data
        dialog = QFileDialog()
        folder = dialog.getExistingDirectory(None, "Select folder to save inside")
        
        if folder is not None:
            # For each energy, save shifted data under unique filename
            if len(sanitized_identifier) > 0:
                filename = f'{sanitized_identifier}_dIdVmap'
            else:
                filename = 'dIdVmap'
            for E, data in self.energy_dict.items():
                # Shift data
                for y, x in self.selected_pixels.items():
                    row = data[y].copy()
                    data[y] = roll(row, -x-1)
                header = f'File format = ASCII\nx-pixels = {data.shape[1]}\ny-pixels = {data.shape[0]}\nEnergy = {E}\nEnergy units = V\ndI/dV maps\n\nStart of data:'
                savetxt(f'{folder}\{filename}_{E}V.asc', data, delimiter='\t', header=header, fmt='%.6f')
            
            # Save point coordinates of slice into '.asc' file
            pts_header = f'File format = ASCII\nStart of data:'
            if len(sanitized_identifier) > 0:
                pts_filename = f'{sanitized_identifier}_glitch_pts.asc'
            else:
                pts_filename = f'glitch_pts.asc'
            pts = fliplr(array(list(self.selected_pixels.items()))) # Reformat dictionary with entries 'y: x' into numpy array with entries [x, y]
            savetxt(f'{folder}\{pts_filename}', pts, header=pts_header, delimiter='\t', fmt='%.6f')
            
            # Show saved message in status bar
            save_message = 'All dI/dV maps'
            save_message += f' saved at {ctime()}'
            self.status_bar.showMessage(save_message, int(1e4))

# Application execution
if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # Show the prompt dialog to select files
    init_prompt = PromptDialog()
    init_prompt.exec_()
    
    # Load and check data
    data = Data3D(init_prompt.file_names)
    if not data.valid:
        print('Data not compatible.')
        sys.exit()
    
    # Load and check file with slice coordinates
    pts = []
    try:
        pts = loadtxt(init_prompt.pt_name, delimiter='\t')
    except:
        pass
    
    # Show main window and execute
    window = MainWindow(data, pts)
    window.show()
    
    sys.exit(app.exec_())