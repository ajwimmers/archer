# ARCHER

Automated Rotational Center Hurricane Eye Retrieval

## Notes

1. This package works with the standard Anaconda Python3 distribution,
	with one exception: it also needs the pyproj module. One way to add 
	this module is

'''bash
conda install pyproj
'''

2. The script demo.py has numerous examples for testing the algorithm
	on a variety of input data. Follow the comments at the beginning 
	of the script to modify this file to initiate tests from various inputs 
	from the data/ subdirectory. By default, demo.py tests nothing, so 
	you'll need to edit it to do anything.

3. The file archer/archer.py contains the *long* header that lays out all 
	the components of input and output. It also contains the rules of 
	proper usage. As user-friendly as I try to make this, there are a lot
	of behaviors of center-fixing (parallax, image format, algorithm
	dependencies on radiation type) that one needs to be in control of 
	in order to get the expected answer.

4. This package is in beta mode. I am open to requests, including:

	a. Better error-trapping and elegant exiting from improper/incomplete 
		data inputs
	b. Better/different status messages during execution, or messages that 
		follow a GeoIPS standard
	c. Capabilities to input other data sources not included (yet) in
		demo.py
