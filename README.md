# fits-movie
Scripts to make a movie pan through of a fits image


Use the path class in 'make_movie.py' to define a path to travel through a fits image, by defining panto's, goto's, dwell's and zoom's. This then calls 'make_field_plot.py' to make the individual frames of the movie. Adjust the 'plot_image' function in 'make_field_plot.py' to make the frames as you wish.

Notes: every frame needs to fit completely in the fits image, it will give errors if this is not the case.
