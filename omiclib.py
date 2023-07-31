import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
from matplotlib.markers import MarkerStyle
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.use("Qt5Agg")

from skimage.measure import profile_line
from qtrangeslider import QRangeSlider
from time import time
import scipy.fftpack as sfft
import numpy as np
import math
import access2thematrix

from PyQt5 import QtWidgets, QtGui, QtCore

'''
Functions in library:
- angle : Returns angle in degrees between +x-axis and points
- find_nearest : Returns the index of the closest element in an array to the given value
- format_file_names : Returns a string condensing an array of file names

Classes in library:
(Data Handlers)
- Data2D : handle 2D data from Scienta-Omicron mtrx files or user-generated .asc files
- Data3D : handle 3D data from Scienta-Omicrom mtrx files or user-generated .asc files

- DraggablePoint : Dynamically updating point that can be dragged by user around graph, in pairs will create a line

(Dynamically updating widgets)
- SliderAndFloatInput : Slider that updates with a float textbox and two buttons to increase and decrease by 1 increment
- RangeSliderAndFloatInputs : Ranged slider with two float textboxes that set the slider range
- PointSetter : 4 float textboxes that set the position of two DraggablePoints

(Canvases)
- LineCutGraph : Line cut performed on topo data, a topo plot with two dynamically updating DraggablePoints and a 1D line cut plot
- dIdVGraph : dI/dV map that can be updated when an updated energy and data is specified
- FFTGraph : FFT of dI/dV map of any given energy and data 
- FFTLineCutGraph : Line cut performed on FFT of dI/dV map and QPI plot
- GlitchFixGraph : Interactive plot where user can select pixels in each row to shift data
'''

''' Standalone functions '''

def angle(x_0, y_0, x_f, y_f):
    # Return the angle in degrees from the +x-axis between pt0 and pt1
    deltaY = y_f - y_0
    deltaX = x_f - x_0
    return math.degrees(math.atan2(deltaY, deltaX))

def find_nearest(array, value):
    # Return the index of the closest element in an array to the given value
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx - 1
    else:
        return idx
    
def format_file_names(file_names):
    # Returns a string for an array of file names, showing the first and last 5 files and the last 25 characters of each
    display_names = [f'...{name[-25:]}' if len(name) > 25 else name for name in file_names]
    first_five = '\n'.join(name.rjust(25, '.') for name in display_names[:5])
    last_five = '\n'.join(name.rjust(25, '.') for name in display_names[-5:])
    if len(file_names) <= 10:
        return '\n'.join(display_names) + f'\n\n{len(file_names)} selected files.'
    return f'{first_five}\n...\n{last_five}\n\n{len(file_names)} selected files.'

''' Scienta-Omicron mtrx data handlers '''

class Data2D:
    def __init__(self, file_name):
        # Initialize the Data2D object by reading the data from the specified file.
        self.valid = True
        self.file_name = file_name
        self.header = ''
        try:
            # Check the file extension to decide how to read the data.
            if file_name[-3:] == 'asc' or file_name[-3:] == 'txt':
                self.txt_to_array()
            elif file_name[-4:] == 'mtrx':
                self.mtrx_to_array()
            else:
                self.valid = False
        except:
            self.valid = False

    def txt_to_array(self):
        # Read data from a text file (ASCII format) and extract header information.
        self.data = np.loadtxt(self.file_name, delimiter='\t')
        with open(self.file_name) as f:
            for line in f:
                if line.startswith('#'):
                    self.header += line
                    # Extract x and y lengths from the header to compute coordinate arrays.
                    if 'x-length' in line:
                        first, remainder = line.split('x-length = ')
                        self.x = float(remainder.split()[0]) / 2
                        self.x = np.linspace(-self.x, self.x, self.data.shape[1])
                    if 'y-length' in line:
                        first, remainder = line.split('y-length = ')
                        self.y = float(remainder.split()[0]) / 2
                        self.y = np.linspace(-self.y, self.y, self.data.shape[0])
        self.header = self.header[:-1]
    
    def mtrx_to_array(self):
        # Read data from a .mtrx file (matrix format) using the access2thematrix library.
        mtrx_data = access2thematrix.MtrxData()
        traces, message = mtrx_data.open(self.file_name)
        selected_image, message = mtrx_data.select_image(traces[0])
        # Convert the data to nanometers and create coordinate arrays.
        self.data = selected_image.data * 1e9
        self.x = selected_image.width * 0.5e9
        self.y = selected_image.height * 0.5e9
        self.x = np.linspace(-self.x, self.x, self.data.shape[1])
        self.y = np.linspace(-self.y, self.y, self.data.shape[0])
        # Create a header for the ASCII format.
        self.header = f'File format = ASCII\nx-pixels = {selected_image.data.shape[1]}\ny-pixels = {selected_image.data.shape[0]}\nx-length = {selected_image.width*1e9}\ny-length = {selected_image.height*1e9}\nz-unit = nm\nStart of data:'

class Data3D:
    def __init__(self, file_names):
        self.valid = True
        self.file_names = file_names
        self.header = ''
        self.energy_dict = {}

        # Check the file extension to determine the file format and read the data accordingly
        if self.file_names[0].endswith('.Aux2(V)_mtrx'):
            self.read_mtrx_files()
        elif self.file_names[0].endswith('.asc'):
            self.read_asc_files()
        else:
            self.valid = False

    def read_mtrx_files(self):
        try:
            # Read data from .Aux2(V)_mtrx files using the access2thematrix module
            mtrx_data = access2thematrix.MtrxData()
            mtrx_data.open(self.file_names[0])
            self.data = mtrx_data.volume_scan['forward/up']['trace']
            self.V = np.round(mtrx_data.scan[0], 2)
            
            # Create a dictionary mapping energy values to 2D data arrays
            self.energy_dict = {self.V[i]: self.data[:, :, i] for i in range(self.data.shape[2])}
        except:
            self.valid = False

    def read_asc_files(self):
        self.V = []
        self.data = []
        for f in self.file_names:
            try:
                data = np.loadtxt(f, delimiter='\t')
                with open(f) as fi:
                    for line in fi:
                        if '# Energy = ' in line:
                            # Extract the energy value from the header
                            V = round(float(line.split('Energy = ')[1].split()[0]), 2)
                            self.V.append(V)
                            self.data.append(data)
                            break 
            except:
                continue

        if len(self.V) > 0:
            # Sort the data and energy values based on energy in ascending order
            sorted_indices = np.argsort(self.V)
            self.V = np.array(self.V)[sorted_indices]
            self.data = np.array(self.data)[sorted_indices]
            
            # Create a dictionary mapping energy values to 2D data arrays
            self.energy_dict = dict(zip(self.V, self.data))
        else:
            # If no valid data was found in the .asc files, set valid to False
            self.valid = False

''' Draggable point class '''

from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
from matplotlib.markers import MarkerStyle

class DraggablePoint:
    # A class for creating draggable points on a Matplotlib plot. The points can be dragged interactively to update a line cut plot.
    # https://stackoverflow.com/questions/28001655/draggable-line-with-draggable-points
    lock = None  # Only one point can be animated at a time
    def __init__(self, parent, x=0.1, y=0.1, size=0.1):
        # Initialize the DraggablePoint object.
        self.parent = parent
        self.point = patches.Ellipse((x, y), size, size, fc='deeppink', alpha=0.65, edgecolor='deeppink', zorder=10)
        self.line = None
        self.x = x
        self.y = y
        parent.fig.axes[0].add_patch(self.point)
        self.press = None
        # self.background = None
        self.connect()

        if self.parent.list_points:
            line_x = [self.parent.list_points[0].x, self.x]
            line_y = [self.parent.list_points[0].y, self.y]
            t = Affine2D().rotate_deg(angle(line_x[0], line_y[0], line_x[1], line_y[1]))
            m = MarkerStyle('>', 'left', transform=t)
            self.line = Line2D(line_x, line_y, color='deeppink', alpha=0.65, zorder=10, marker=m, markevery=(1, 2))
            parent.fig.axes[0].add_line(self.line)

    def connect(self):
        # Connect events
        self.cidpress = self.point.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.point.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion = self.point.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        # Handle the button press event for the draggable point.
        if event.inaxes != self.point.axes:
            return
        if DraggablePoint.lock is not None:
            return
        contains, attrd = self.point.contains(event)
        if not contains:
            return
        self.press = (self.point.center), event.xdata, event.ydata
        DraggablePoint.lock = self

        # Draw everything but the selected point and store the pixel buffer
        canvas = self.point.figure.canvas
        axes = self.point.axes
        self.point.set_animated(True)
        if self == self.parent.list_points[1]:
            self.line.set_animated(True)
        else:
            self.parent.list_points[1].line.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.point.axes.bbox)

        # Redraw just the point
        axes.draw_artist(self.point)

        # Blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        # Handle the mouse motion event when dragging the point.
        if DraggablePoint.lock is not self:
            return
        if event.inaxes != self.point.axes:
            return
        self.point.center, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress

        self.point.center = (self.point.center[0] + dx, self.point.center[1] + dy)

        # Update the line cut plot and redraw
        self.parent.updateLineCut()
        self.parent.draw()

        canvas = self.point.figure.canvas
        axes = self.point.axes
        # Restore the background region
        canvas.restore_region(self.background)

        # Redraw the point and line
        axes.draw_artist(self.point)
        if self == self.parent.list_points[1]:
            axes.draw_artist(self.line)
        else:
            axes.draw_artist(self.parent.list_points[1].line)

        # Blit just the redrawn area
        canvas.blit(axes.bbox)

        self.x = self.point.center[0]
        self.y = self.point.center[1]
        try:
            self.parent.parent.updateCoords()
        except:
            pass

        if self == self.parent.list_points[1]:
            line_x = [self.parent.list_points[0].x, self.x]
            line_y = [self.parent.list_points[0].y, self.y]
            t = Affine2D().rotate_deg(angle(line_x[0], line_y[0], line_x[1], line_y[1])) # Rotation to the angle of the line
            self.line.set_marker(MarkerStyle('>', 'left', transform=t)) # Marker to turn the line into an arrow
            self.line.set_data(line_x, line_y)
        else:
            line_x = [self.x, self.parent.list_points[1].x]
            line_y = [self.y, self.parent.list_points[1].y]
            t = Affine2D().rotate_deg(angle(line_x[0], line_y[0], line_x[1], line_y[1]))
            self.parent.list_points[1].line.set_marker(MarkerStyle('>', 'left', transform=t))
            self.parent.list_points[1].line.set_data(line_x, line_y)

        # Blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        # Handle the button release event when dragging the point.
        if DraggablePoint.lock is not self:
            return

        self.press = None
        DraggablePoint.lock = None

        # Turn off the animation property and reset the background
        self.point.set_animated(False)
        if self == self.parent.list_points[1]:
            self.line.set_animated(False)
        else:
            self.parent.list_points[1].line.set_animated(False)

        self.background = None

        # Redraw the full figure
        self.point.figure.canvas.draw()

        self.x = self.point.center[0]
        self.y = self.point.center[1]

    def disconnect(self):
        # Disconnect all the stored connection ids.
        self.point.figure.canvas.mpl_disconnect(self.cidpress)
        self.point.figure.canvas.mpl_disconnect(self.cidrelease)
        self.point.figure.canvas.mpl_disconnect(self.cidmotion)

    def delete(self):
        # Remove DraggablePoint from the canvas.
        if self.parent is None:
            return
        self.disconnect()
        self.point.remove()
        if self.line is not None:
            self.line.remove()

        # Redraw the figure
        self.parent.draw_idle()
        self.parent = None

''' Dynamically updating widgets '''
            
# Used for energy sliders, creates a slider and a textbox which takes float inputs that update dynamically
class SliderAndFloatInput(QtWidgets.QWidget):
    def __init__(self, title, V, parent=None):
        super().__init__()
        
        self.parent = parent
        self.V = np.sort(V)
        self.V_step = np.round(np.abs(V[1] - V[0]), 2)
        self.V_min = V[0]
        self.V_max = V[-1]
        
        # Create slider, textbox, and increase and decrease buttons
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setRange(0, V.size - 1)  # Set the range of the slider, can only take integer values
        self.slider.setSingleStep(1)  # Set the step size for slider movements
        self.slider.setPageStep(1)
        
        self.float_textbox = QtWidgets.QLineEdit()
        self.float_textbox.setValidator(QtGui.QDoubleValidator())  # Allow only float values
        self.float_textbox.setText('{:g}'.format(float('{:.2g}'.format(V[0]))))  # Set the initial value of the input text
        
        self.decrease_button = QtWidgets.QPushButton("<<")
        self.increase_button = QtWidgets.QPushButton(">>")
        
        # Connect signals and slots
        self.float_textbox.textEdited.connect(self.updateSlider)
        self.slider.valueChanged.connect(self.updateFloatTextbox)
        self.decrease_button.clicked.connect(self.decreaseValue)
        self.increase_button.clicked.connect(self.increaseValue)
        
        # Arrange widget in layouts
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.decrease_button)
        layout.addWidget(self.float_textbox)
        layout.addWidget(self.increase_button)
        main_layout = QtWidgets.QVBoxLayout(self)
        main_layout.addWidget(QtWidgets.QLabel(title))
        main_layout.addLayout(layout)
        main_layout.addWidget(self.slider)
        
        self.setLayout(main_layout)

    def updateFloatTextbox(self, value):
        float_value = self.V[value]
        self.float_textbox.setText('{:g}'.format(float('{:.2g}'.format(float_value))))
        self.parent.updatePlots(float_value)
        
    def updateSlider(self, text):
        try:
            float_value = round(float(text), 2)
            float_value = max(self.V_min, min(self.V_max, float_value))  # Clamp the value between V_min and V_max
            slider_value = np.argmin(np.abs(self.V - float_value)) # Find index of V closest to the given value
            self.slider.setValue(slider_value)
        except ValueError:
            pass
        
    def decreaseValue(self):
        current_value = self.float_textbox.text()
        try:
            float_value = float(current_value)
            float_value -= self.V_step
            float_value = round(max(self.V_min, float_value), 2)  # Clamp the value to V_min
            self.float_textbox.setText('{:g}'.format(float('{:.2g}'.format(float_value))))
            self.updateSlider(str(float_value))
        except ValueError:
            pass
        
    def increaseValue(self):
        current_value = self.float_textbox.text()
        try:
            float_value = float(current_value)
            float_value += self.V_step
            float_value = round(min(self.V_max, float_value), 2)  # Clamp the value to V_max
            self.float_textbox.setText('{:g}'.format(float('{:.2g}'.format(float_value))))
            self.updateSlider(str(float_value))
        except ValueError:
            pass

# Used for adjusting colormap, range slider with float inputs for max and min values 
class RangeSliderAndFloatInputs(QtWidgets.QWidget):
    def __init__(self, title, min_val, max_val, parent=None, slider_id=1):
        super().__init__(parent)
        self.main_window = parent
        self.id = slider_id
        
        # Create a linearly spaced array of values between min_val and max_val for mapping slider positions to float values
        self.map = np.round(np.linspace(min_val, max_val, 100), 2)
        
        # Create a range slider widget
        self.range_slider = QRangeSlider(QtCore.Qt.Horizontal)
        self.range_slider.setRange(0, 99)
        self.range_slider.setValue((0, 49))  # Set the initial range of the slider
        
        # Create text boxes to display the float values corresponding to the slider positions
        self.left_textbox = QtWidgets.QLineEdit(str(self.map[0]))
        self.right_textbox = QtWidgets.QLineEdit(str(self.map[49]))
        
        # Connect signals and slots for handling value changes
        self.range_slider.valueChanged.connect(self.updateTextboxes)
        self.left_textbox.setValidator(QtGui.QDoubleValidator()) 
        self.right_textbox.setValidator(QtGui.QDoubleValidator())
        self.left_textbox.editingFinished.connect(self.updateRangeSliderFromTextbox)
        self.right_textbox.editingFinished.connect(self.updateRangeSliderFromTextbox)
        
        # Create the layout for the widget
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(QtWidgets.QLabel(title))  
        layout.addWidget(self.range_slider) 
        textbox_layout = QtWidgets.QHBoxLayout()
        textbox_layout.addWidget(self.left_textbox)
        textbox_layout.addWidget(self.right_textbox)
        layout.addLayout(textbox_layout)
        
        self.setLayout(layout)
        
    def updateTextboxes(self, vals):
        # Update the text boxes with the float values corresponding to the slider positions
        lower = self.map[vals[0]]
        upper = self.map[vals[1]]
        self.left_textbox.setText(str(lower))
        self.right_textbox.setText(str(upper))
        
        # Notify the main window about the range change
        self.main_window.updateCmap(lower, upper, self.id)
        
    def updateRangeSliderFromTextbox(self):
        # Update the range slider positions based on the values in the text boxes
        lower = float(self.left_textbox.text())
        upper = float(self.right_textbox.text())
        
        # Find the closest values in the map array to the entered float values
        lower_index = np.searchsorted(self.map, lower)
        upper_index = np.searchsorted(self.map, upper)
        
        # Set the slider positions based on the found indices
        self.range_slider.setValue((lower_index, upper_index))
        
        # Notify the main window about the range change
        self.main_window.updateCmap(self.map[lower_index], self.map[upper_index], self.id)
        
    def getRange(self):
        # Return the current range as floats
        return float(self.left_textbox.text()), float(self.right_textbox.text())
    
class PointSetter(QtWidgets.QWidget):
    def __init__(self, graph):
        super().__init__()
        self.graph = graph
        self.labels = ['Start X:', 'Start Y:', 'End X:', 'End Y:']
        self.textboxes = [QtWidgets.QLineEdit() for _ in range(4)]

        # Set validators and placeholder text for each textbox
        for i, tb in enumerate(self.textboxes):
            tb.setValidator(QtGui.QDoubleValidator())
            tb.setPlaceholderText(self.labels[i].split(':')[0])

        # Use a grid layout for arrangement of labels and textboxes
        layout = QtWidgets.QGridLayout()
        for row in range(2):
            for col in range(2):
                label = QtWidgets.QLabel(self.labels[row * 2 + col])
                layout.addWidget(label, row, col * 2)
                textbox = self.textboxes[row * 2 + col]
                textbox.setObjectName(f'textbox_{self.labels[row * 2 + col]}')
                layout.addWidget(textbox, row, col * 2 + 1)

        self.setLayout(layout)

        # Connect editingFinished signals to update functions
        self.textboxes[0].editingFinished.connect(self.updatePoint)
        self.textboxes[1].editingFinished.connect(self.updatePoint)
        self.textboxes[2].editingFinished.connect(self.updatePoint)
        self.textboxes[3].editingFinished.connect(self.updatePoint)

    def updatePoint(self):
        values = []
        for textbox in self.textboxes:
            text = textbox.text()
            if not text:
                values.append(None) # If a textbox is empty, use the corresponding point's center coordinates
            else:
                try:
                    values.append(float(text))
                except ValueError:
                    return

        # Check if any textbox is empty and update the values accordingly
        pt1 = (values[0] if values[0] is not None else self.graph.list_points[0].point.center[0],
            values[1] if values[1] is not None else self.graph.list_points[0].point.center[1])

        pt2 = (values[2] if values[2] is not None else self.graph.list_points[1].point.center[0],
            values[3] if values[3] is not None else self.graph.list_points[1].point.center[1])

        self.graph.updatePlots(pt1, pt2)
        self.graph.draw_idle()

''' Graphs '''

       
class LineCutGraph(FigureCanvas):
    # Topo line cuts, one plot with topography and DraggablePoints as well as line cut plot
    def __init__(self, parent=None, width=4, height=5, dpi=150):
        super().__init__()
        self.parent = parent
        # Initialize the figure and subplot
        self.fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.fig.subplots_adjust(hspace=0.5)
        self.topo_ax = self.fig.add_subplot(3, 1, (1, 2), aspect='equal')
        self.topo_ax.set_title('STM Topography', fontsize=8)
        self.topo_ax.set_xlabel('x (nm)', fontsize=7)
        self.topo_ax.set_ylabel('y (nm)', fontsize=7)
        self.topo_ax.tick_params(axis='both', labelsize=7)
        self.lc_ax = self.fig.add_subplot(313)
        self.lc_ax.set_title('Line cut', fontsize=8)
        self.lc_ax.set_xlabel('x (nm)', fontsize=7)
        self.lc_ax.set_ylabel('z (nm)', fontsize=7)
        self.lc_ax.tick_params(axis='both', labelsize=7)
        self.line = Line2D([], [], c='deeppink', lw=1)
        self.lc_ax.add_line(self.line)
        self.im = None
        self.c = self.fig.colorbar(self.im, ax=self.topo_ax, label='z (nm)', fraction=0.02)
        self.c.ax.tick_params(labelsize=7)
        self.data = None
        self.x = None
        self.y = None
        self.header = None
        self.list_points = []

        # Initialize the FigureCanvas
        FigureCanvas.__init__(self, self.fig)

        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def updatePlots(self, pt1, pt2):
        # Update the STM Topography and line cut plots
        self.im = self.topo_ax.imshow(self.data, cmap='Blues_r', origin='lower',
                                      extent=[self.x[0], self.x[-1], self.y[0], self.y[-1]], zorder=1)
        self.c.mappable.set_clim(self.im.get_clim())
        self.c.mappable.set_cmap('Blues_r')
        self.plotDraggablePoints(pt1, pt2, self.x[-1] / 20)
        self.draw_idle()

    def plotDraggablePoints(self, xy1, xy2, size=0.1):
        # Hide draggable points that are already shown
        for l in self.list_points:
            l.delete()
        self.list_points = []

        # Add new draggable points at desired coordinates
        self.list_points.append(DraggablePoint(self, xy1[0], xy1[1], size))
        self.list_points.append(DraggablePoint(self, xy2[0], xy2[1], size))

        # Update line cut and QPI plots based on new line cut
        self.updateLineCut()

    def updateLineCut(self):
        # Update the line cut based on the draggable points' position
        if self.data is None:
            return

        pt1_xy = self.list_points[0].point.center
        pt2_xy = self.list_points[1].point.center
        # Find the indices of the x and y coordinates
        pt1 = [find_nearest(self.y, self.list_points[0].point.center[1]),
               find_nearest(self.x, self.list_points[0].point.center[0])]
        pt2 = [find_nearest(self.y, self.list_points[1].point.center[1]),
               find_nearest(self.x, self.list_points[1].point.center[0])]

        line_cut_data = profile_line(self.data, pt1, pt2, order=1, mode='constant')  # Returns interpolated data along given line
        # Plot and resize
        dist = math.dist(pt1_xy, pt2_xy)
        self.line.set_data(np.linspace(0, dist, line_cut_data.size), line_cut_data)
        self.lc_ax.set_xlim(0, dist)
        self.lc_ax.set_ylim(min(line_cut_data), max(line_cut_data))

    def savePlot(self):
        # Save the LineCutGraph plot to a file
        file_path, _ = QtWidgets.QFileDialog.getSaveFileName(None, 'Save Plot', '', 'EPS Files (*.eps);;PNG Files (*.png);;All Files (*)')
        if file_path:
            self.fig.savefig(file_path, dpi=300, facecolor='w')
            file_path = file_path[:-3] + 'asc'
            header = f'File format = ASCII\nStart point = {self.list_points[0].point.center}\nEnd point = {self.list_points[1].point.center}\n\n\nz-unit = nm\nStart of data:'
            np.savetxt(file_path, list(zip(*self.line.get_data())), delimiter='\t', header=header, fmt='%.6f')
            print('Plot saved!')
        else:
            print('No file path selected.')

class dIdVGraph(FigureCanvas):
    # colormap of dI/dV map given data and energy
    def __init__(self, data, energy, parent=None, dpi=150):
        super().__init__()

        self.fig = Figure(figsize=(4, 4), dpi=dpi, tight_layout=True)        
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_title(f'dI/dV Map {energy}V')
        self.im = self.ax.imshow(data, origin='lower', cmap='Blues_r')
        self.im.set_clim(np.mean(data) - 2 * np.std(data), np.mean(data) + 2 * np.std(data))
        self.ax.set_xlabel('x (pixels)')
        self.ax.set_ylabel('y (pixels)')

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def updateFigure(self, data, energy):
        # Update the image data using set_data
        self.im.set_data(data)
        self.ax.set_title(f'dI/dV Map {energy}V')
        self.im.set_clim(np.mean(data) - 2 * np.std(data), np.mean(data) + 2 * np.std(data))
        self.draw_idle()
   
class FFTGraph(FigureCanvas):
    # colormap of FFT of dI/dV map
    def __init__(self, data, parent=None, dpi=150):
        super().__init__()
        
        self.fig = Figure(figsize=(4, 4), dpi=dpi, tight_layout=True)
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_title('FFT')
        self.ax.set_xlabel('$k_x$')
        self.ax.set_ylabel('$k_y$')
        
        # Take 2D-FFT and center around the origin
        data_fft = sfft.fft2(data)
        self.data_shift = np.abs(sfft.fftshift(data_fft))
        x_fft = sfft.fftshift(np.fft.fftfreq(data.shape[1], d=1))
        y_fft = sfft.fftshift(np.fft.fftfreq(data.shape[0], d=1))
        self.im = self.ax.imshow(self.data_shift, origin='lower', cmap='bone_r', vmin=0, vmax=np.mean(self.data_shift) + 2 * np.std(self.data_shift), extent=[x_fft[0], x_fft[-1], y_fft[0], y_fft[-1]])        

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def getMax(self):
        # Return a reasonable maximum value for the colormap
        return np.mean(self.data_shift) + 4 * np.std(self.data_shift)
    
    def updateFigure(self, data):
        # Update plot when data is changed
        data_fft = sfft.fft2(data)
        self.data_shift = np.abs(sfft.fftshift(data_fft))
        self.im.set_data(self.data_shift)
        self.draw_idle()
        
    def updateClim(self, lower, upper):
        # Update limits of colormap on the plot
        self.im.set_clim(vmin=lower, vmax=upper)
        self.draw_idle()
 
class FFTLineCutGraph(FigureCanvas):
    # Canvas with FFT, line cut, and QPI plots
    def __init__(self, data, current_energy, parent=None, dpi=150):
        super().__init__()
        self.parent = parent
        
        # Take 2D-FFT of all energies and store in dictionary
        self.fft_dict = {data.V[i]: np.abs(sfft.fftshift(sfft.fft2(data.data[:,:,i]))) for i in range(len(data.V))}
        self.x_fft = sfft.fftshift(np.fft.fftfreq(data.data.shape[1], d=1))
        self.y_fft = sfft.fftshift(np.fft.fftfreq(data.data.shape[0], d=1))
        self.current_energy = current_energy
        
        # Initialize axes
        space = 0.07
        self.fig = Figure(figsize=(4, 8), dpi=dpi)
        self.fft_ax = self.fig.add_axes([0 + space, 0.4 + space, 0.6 - 2 * space, 0.6 - 2 * space], aspect='equal')
        self.lc_ax = self.fig.add_axes([0 + 1.25 * space, 0 + 1.25 * space, 0.6 - 2.5 * space, 0.4 - 2.5 * space])
        self.qpi_ax = self.fig.add_axes([0.6 + .25 * space, 0 + .75 * space, 0.4 - 1.5 * space, 1 - 1.5 * space], aspect='equal')
        self.fig.subplots_adjust(hspace=0.5)
        
        # Initialize FFT and line cut plots
        self.fft_ax.set_title('FFT', fontsize=8)
        self.fft_ax.set_xlabel('$k_x$', fontsize=7, labelpad=0)
        self.fft_ax.set_ylabel('$k_y$', fontsize=7, labelpad=0)
        self.fft_ax.tick_params(axis='both', labelsize=7)
        
        self.lc_ax.set_title('Line cut', fontsize=8)
        self.lc_ax.set_xlabel('$k_x$', fontsize=7, labelpad=0)
        self.lc_ax.set_ylabel('FFT', fontsize=7, labelpad=0)
        self.lc_ax.tick_params(axis='both', labelsize=7)
        
        self.fft_im = self.fft_ax.imshow(self.fft_dict.get(self.current_energy), origin='lower', cmap='bone_r', extent=[self.x_fft[0], self.x_fft[-1], self.y_fft[0], self.y_fft[-1]], vmin=0, vmax=self.getMax() / 2)
        self.line = Line2D([], [], c='deeppink', lw=1)
        self.lc_ax.add_line(self.line)
        self.lc_ax.set_ylim(self.fft_im.get_clim()[0], self.fft_im.get_clim()[1])
        
        # Store draggable points and initialize the QPI plot
        self.list_points = []
        self.updatePlots([self.x_fft.max() * 2 / 3, self.y_fft.max() * 2 / 3],
                                        [self.x_fft.min() * 2 / 3, self.y_fft.min() * 2 / 3])
        self.c = self.fig.colorbar(self.qpi_im, ax=self.qpi_ax, label='FFT', fraction=0.025)
        self.c.ax.tick_params(labelsize=7)

        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        self.draw_idle()

    def updatePlots(self, pt1, pt2):
        # Resets QPI plot which has to be redrawn, and calls a function that will plot new data
        self.qpi_ax.clear()
        data = self.fft_dict.get(self.current_energy)
        self.fft_im.set_data(data)
        
        self.qpi_ax.set_title('QPI', fontsize=8)
        self.qpi_ax.set_xlabel(r'$\vec{q}$', fontsize=7, labelpad=0)
        self.qpi_ax.set_ylabel('E', fontsize=7, labelpad=0)
        self.qpi_ax.tick_params(axis='both', labelsize=7)

        self.plotDraggablePoints(pt1, pt2, self.x_fft.max() / 20)
        
    def plotDraggablePoints(self, xy1, xy2, size=0.1):
        # Hide draggable points that are already shown
        for l in self.list_points:
            l.delete()
        self.list_points = []
        
        # Add new draggable points at desired coordinates
        self.list_points.append(DraggablePoint(self, xy1[0], xy1[1], size))
        self.list_points.append(DraggablePoint(self, xy2[0], xy2[1], size))
        
        # Update line cut and QPI plots based on new line cut
        self.updateLineCut()

    def updateLineCut(self):
        # Find points in FFT closest to given points
        pt1_xy = self.list_points[0].point.center
        pt2_xy = self.list_points[1].point.center
        pt1 = [find_nearest(self.y_fft, self.list_points[0].point.center[1]),
               find_nearest(self.x_fft, self.list_points[0].point.center[0])]
        pt2 = [find_nearest(self.y_fft, self.list_points[1].point.center[1]),
               find_nearest(self.x_fft, self.list_points[1].point.center[0])]
        
        # Find the linear interpolation along line for all energies
        self.all_lines = {}
        for en, dat in self.fft_dict.items():
            self.all_lines[en] = profile_line(dat, pt1, pt2, order=1, mode='constant')

        # Update line cut plot
        current_line = self.all_lines.get(self.current_energy)
        dist = math.dist(pt1_xy, pt2_xy)
        x = np.linspace(0, dist, current_line.size)
        self.line.set_data(x, current_line)
        self.lc_ax.set_xlim(0, dist)
        
        # Update QPI plot
        E, self.z = zip(*sorted(self.all_lines.items()))
        self.qpi_ax.set_xlim(0, dist)
        self.qpi_im = self.qpi_ax.pcolormesh(x, E, self.z, vmin=0, vmax=self.getMax() / 2)
        
    def updateClimFFTLC(self, lower, upper):
        # Update colormap limits for FFT and line cut plots
        self.fft_im.set_clim(vmin=lower, vmax=upper)
        self.lc_ax.set_ylim(lower, upper)
        
    def updateClimQPI(self, lower, upper):
        # Update colormap limits for QPI plot
        self.qpi_im.set_clim(vmin=lower, vmax=upper)
        self.c.mappable.set_clim(vmin=lower, vmax=upper)
        
    def updateFigure(self, current_energy):
        # If the energy is changed, updates all plots accordingly
        vmin2, vmax2 = self.qpi_im.get_clim()
        self.current_energy = current_energy
        self.updatePlots(self.list_points[0].point.center, self.list_points[1].point.center)
        self.updateClimQPI(vmin2, vmax2)
        self.draw_idle()
    
    def getMax(self):
        # Returns reasonable maximum 
        data = list(self.fft_dict.values())
        return np.mean(data) + 5 * np.std(data)
    
    def saveData(self, folder, identifier=''):
        # Save line cut and QPI data
        lc = np.array(self.line.get_data(orig=False)).transpose() # Returns pairs of points
        qpi = self.z
        
        # Set headers and filenames
        E = self.current_energy
        lc_header = f'File format = ASCII\nStart of data:'
        qpi_header = f'File format = ASCII\nx-pixels = {lc.shape[0]}\ny-pixels = {len(self.fft_dict)}\nx-length = {lc[:,0].max()}\nmin E (y) = {min(self.fft_dict)}\nmax E (y) = {max(self.fft_dict)}\nE (y) units = V\nStart of data:'
        if len(identifier) > 0:
            lc_filename = f'{identifier}_line_cut_{E}V.asc'
            qpi_filename = f'{identifier}_qpi.asc'
        else:
            lc_filename = f'line_cut_{E}V.asc'
            qpi_filename = 'qpi.asc'
            
        # Save as '.asc' files
        np.savetxt(f'{folder}\{lc_filename}', lc, delimiter='\t', header=lc_header, fmt='%.6f')
        np.savetxt(f'{folder}\{qpi_filename}', qpi, delimiter='\t', header=qpi_header, fmt='%.6f')

class GlitchFixGraph(FigureCanvas):
    # Interactive canvas to select pixels and shift data
    def __init__(self, pts, parent=None):
        # Initialize the figure and the canvas
        self.fig = Figure()
        super().__init__(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.parent = parent
        self.mpl_connect('button_press_event', self.on_double_click)

        # Initialize the data and original_data attributes
        self.data = self.parent.energy_dict[self.parent.current_energy]
        self.original_data = self.data
        self.im = self.ax.imshow(self.data, origin='lower', cmap='Blues_r')  # Display the data as an image
        self.ax.set_xlabel('x (pixels)')
        self.ax.set_ylabel('y (pixels)')
        self.ax.set_title(f'dI/dV Map {self.parent.current_energy}V')

        # Initialize the selected_pixels attribute
        self.selected_pixels = {}
        if len(pts) > 0:
            self.selected_pixels = {int(x[1]): int(x[0]) for x in pts}

        # Initialize red_circles to store circle patches representing selected pixels
        self.red_circles = {}
        for y, x in self.selected_pixels.items():
            circle = plt.Circle((x, y), 0.25, fc='red', alpha=0.75, edgecolor='red', zorder=10)
            self.red_circles[y] = circle
            self.ax.add_patch(circle)  # Add the circle patch to the plot

        self.ax.set_xlim(0, self.data.shape[1])
        self.ax.set_ylim(0, self.data.shape[0])

    def update_data(self, current_energy):
        # Update the data and original_data attributes based on the current_energy
        self.data = self.parent.energy_dict[current_energy]
        self.original_data = self.parent.energy_dict[current_energy]
        self.update_plot()

    def update_plot(self):
        # Update the plot based on the current mode (Edit or Shift)
        self.ax.set_title(f'dI/dV Map {self.parent.current_energy}V')
        if self.parent.mode == 'Edit':
            # Update the image data and color limits in Edit mode
            self.im.set_data(self.original_data)
            self.im.set_clim(np.mean(self.original_data) - 2 * np.std(self.original_data), np.mean(self.original_data) + 2 * np.std(self.original_data))

            # Remove red circles for deselected pixels
            to_remove = set(self.red_circles.keys()).difference(self.selected_pixels)
            for y in to_remove:
                circle = self.red_circles.pop(y)
                circle.remove()  # Remove the circle patch from the plot

            # Add or update red circles for selected pixels
            for y, x in self.selected_pixels.items():
                circle = self.red_circles.get(y)
                if circle is None:
                    circle = plt.Circle((x, y), 0.25, fc='red', alpha=0.75, edgecolor='red', zorder=10)
                    self.red_circles[y] = circle
                    self.ax.add_patch(circle)  # Add the circle patch to the plot
                else:
                    circle.center = (x, y)  # Update the circle's center position

        else:
            # Update the image data and color limits in Shift mode
            self.im.set_data(self.data)
            self.im.set_clim(np.mean(self.data) - 2 * np.std(self.data), np.mean(self.data) + 2 * np.std(self.data))

            # Clear all red circles in Shift mode
            self.clear_red_circles()

        # Update the canvas
        self.fig.canvas.draw_idle()

    def on_double_click(self, event):
        # Handle the double-click event
        if not self.parent.accepts_double_click:
            return
        if event.inaxes != self.ax:
            return
        if event.dblclick:
            x, y = int(event.xdata), int(event.ydata)
            if self.parent.mode == 'Edit':
                # Toggle selection of a pixel in Edit mode
                if y in self.selected_pixels and self.selected_pixels[y] == x:
                    del self.selected_pixels[y]
                else:
                    self.selected_pixels[y] = x
                self.update_plot()
                
    def shift_data(self):
        # Shift the selected pixels' data
        self.fill_selected_pixels()
        self.data = self.original_data.copy()
        for y, x in self.selected_pixels.items():
            row = self.data[y].copy()
            self.data[y] = np.roll(row, -x-1)
            
    def fill_selected_pixels(self):
        # Fill the selected pixels' data based on pixels with lower y-value
        if not self.selected_pixels:
            return
        min_y = min(self.selected_pixels.keys())
        for y in range(min_y + 1, self.data.shape[0]):
            if y in self.selected_pixels:
                continue
            if y - 1 in self.selected_pixels:
                self.selected_pixels[y] = self.selected_pixels[y - 1]
    
    def clear_red_circles(self):
        # Remove all red circles from the plot
        for c in self.red_circles.values():
            c.remove()
        self.red_circles = {}
