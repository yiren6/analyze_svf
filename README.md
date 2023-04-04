# analyze_svf
analyze_svf is an automated program developed in MATLAB for analyzing sooting images and extracting soot volume fraction (SVF) quantitatively. The program can analyze images with both fixed and free-floating droplets during a burn.

## Installation
1. Install MATLAB on your computer.
2. Clone or download the CornellSVF repository from Github.
3. Open MATLAB and navigate to the folder where the repository and image data is located. 
4. Run the program by typing CornellSVF in the MATLAB command window.

## Usage
1. Load the original image and the attenuated image when prompted. The program will automatically convert the images to grayscale intensity matrices with an 8-bit unsigned integer.
2. Measure the droplet radius manually and input a potential range of droplet radius, usually at a range of ±10 pixels from the manual measurement, to facilitate the program to identify the centers of droplets and the radiuses of the droplets for both images.
3. Verify the droplet selection by checking the prompt windows for the images, which will show the droplet boundary marked in a red circle.
4. Input the number of ROI under investigation. For one ROI analysis, the program will ask the user to input the endpoint of analysis using a mouse. For multi-ROI analysis, the user will input a direction-of-analysis (DOA) angle θ by doing two mouse inputs, start and end accordingly.
5. With the DOA angle defined, the program draws other ROIs based on evenly dividing up the DOA angle so that the angle between any two ROIs on the same quadrant is θ/(nROI-1).
6. The program will extract the soot volume fraction for each ROI and output the results.

## Dependencies
CornellSVF uses the following MATLAB toolboxes and functions:
Image Processing Toolbox: imfindcircles.m and improfile.m
Financial Toolbox: tsmovavg.m

## References
If using CornellSVF in your research, please cite the following paper:

[Insert reference to the paper associated with the program here]

## License
MIT license
