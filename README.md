StripesAnalysisSpots is an image analysis platform developed to spatially organize smFish analysis results of drosophile embryos. This software is written in Python™ 3.10.
It is ment to be a post-processing tool of smFish_Software (https://github.com/ant-trullo/smFiSH_software), so it works on the results of analysis obtained with it.

StripesAnalysisSpots was developed on the operating system Linux Ubuntu 24.04 64-bit and tested on Linux Ubuntu 24.04.

Requirements: StripesAnalysisSpots is written in Python3.10, but should work as well with newer versions of Python3. All the packages used by this software (with their version) are listed in the file requirements.txt: to install all of them should be sufficient to run the command 'pip install -r requirements.txt' in your console.

Install  and Run StripesAnalysisSpots: Clone the repository and put all the files in a folder; open with a text editor the file StripesAnalysisSpots.py and modify lines 28, 29, 30, 31, 32 pasting respectively the path of the analysis folder, path to raw data file, specifying 'hor' if you want strips to be organized horizontally or 'vert' if vertically, the tag 'tassels' or 'membrane' if you want to work on pseudo-cells or membrane and finally if the flag '_a' to work with spots_a or '_b' to work with spots_b. Than open a terminal (or cmd for Windows) move in the software folder and run 'python3 StripesAnalysisSpots.py' and press enter. 

A widget will popup and let you manually choose the twwo groups you want.

For any question or issue send an email at: antonio.trullo@igmm.cnrs.fr
