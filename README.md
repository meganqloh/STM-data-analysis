# STM Data Analysis Suite
- `topo-fft.ipynb` shows and saves the FFT of topography data after subtracting a linear background with an adjustable colormap.
- `topo-line-cut.py` allows the user to select line cuts dynamically from topography data by dragging the two endpoints of a line.
- `dI-dV.py` offers three analysis tools for any number of dI/dV maps (with data of the same shape): an interface to switch between all dI/dV maps, dI/dV maps and their FFTs, and the user-selected line cut of FFTs and the QPI for the line cut.
- `glitch-fix.py` uses user input to correct the glitch in dI/dV map data where part of the array is offset.
- `omiclib.py` stores dependencies of `topo-line-cut.py`, `glitch-fix.py`, and `dI-dV.py`. Should be in the same folder as the python files
Created on Python version 3.10.5

## Dependencies
`matplotlib`,
`math`,
`numpy`,
`skimage`,
`time`,
`sys`,
`os`,
`PyQt5`,
`qtrangeslider`,
`pathvalidate`,
`access2thematrix`

(`topo-fft.ipynb` only)
`tkinter`,
`ipywidgets`,
`pyransac3d`,
`pandas`,
`scipy`,
`functools`,
`ipyvuetify`,
`traitlets`,
`re`,
`IPython`

## Acceptable data formats
`.mtrx` files
`.asc` files - for 1D data, as coordinates; for 2D data, as a 2D array. Lines that are not data should start with '#' and data should be delimited with tabs.

## How to run
### `topo-fft.ipynb`
Can be run using Jupyter Notebook or an IDE that supports Jupyter Notebook file formats. `ctrl`+`enter` to run a cell, `shift`+`enter` to run a cell and select the next cell.

### `.py` files
For most systems, the command should be `python3 scriptname.py`. However, running on different systems may have different commands to run a script. The scripts can also be run on any IDEs compatible with Python. Opening the file directly from file explorer in Windows may also work.

### `topo-line-cut.py`
First select the 'Select a file to analyze' button at the top left of the window to upload a `.Z_mtrx` or `.asc` file. Once the data is loaded, you can click and drag the circular endpoints of the line to update the line cut below. The arrow shows the direction of the line cut. The 'Save' button will save an `.eps` file of the screen and an `.asc` file of the line cut as xy-coordinates. The toolbar with symbols can be used to zoom in on the topography plot or change the axes/colormap. The home symbol will reset the plot view. The coordinates of the start and end point of the line cut can also be specified (press enter or unselect the textbox to submit new coordinate).

### `dI-dV.py`
The initialization interface prompts the user to upload file(s) to analyze and choose the function to use. The file(s) to analyze should either be 1 `.Aux2(V)_mtrx` file, or any number of `.asc` files.

In the main window, the file identifier is the unique description to save files under, which will be at the front of every file saved. The `Save all dI/dV maps` button will save all dI/dV maps as `.asc` files with their respective energies in the folder that the user selects. The `Save line cut / QPI` will save the data from the line cut and from the QPI.

### `glitch-fix.py`
The initialization interface has the same upload files to analyze prompt, as well as an optional upload prompt for a file with coordinates to slice along.

In the main window, in `Edit` mode, the user can zoom in with the Matplotlib toolbar and double click a maximum of one pixel per row. The pixel that is selected and every pixel to the left of it will be shifted to the end of the array in `Shift` mode. If no pixel is selected in a row, it will shift the same amount of pixels as the row below it. If a pixel is already selected in a row and the user double clicks another pixel in the row, the original pixel will be unselected and the new pixel will be selected. Pixels can also be unselected by double clicking the selected pixel. The new data will be saved as a group of `.asc` files in the folder that the user selects.

## Future Implementations
- Convert `topo-fft.ipynb` into a `.py` file
- Streamline animation of line cuts

If there are any questions / clarifications, please email meganloh@stanford.edu.
