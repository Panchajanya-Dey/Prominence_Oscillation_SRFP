# These codes are developed by Panchajanya Dey, Indian Institute of Science, Bangalore
# as a part of the SRFP project prominence oscillation in solar corona funded by Indian Academy of Science (2023)
# this folder contains the codes are for downloading the necessary fits files, modifying them by multi-gaussian image processing
# and creating an animation that helps to identify the location of the prominence oscillation
# the next part is to determine the pixel coordinates and creating a 4-pixel wide slit and get the corresponding XT map
# which is also covered in this project
# finally handpicking the prominence oscillation from the XT map and fitting it to a damped harmonic oscillation is done
# the final result is the period, damping time and amplitude determination for a LAL prominence oscillation
# I thank Prof. Vaibhav Pant (ARIES) and his group member Upasna Baweja (ARIES) for helping me throughout the project
# and clearing my doubts
# for comments, suggestions and clarity regarding the code contact at the following email address
# panchajanyad@iisc.ac.in

import sunpy.map
from sunpy.instr.aia import aiaprep
from sunpy.net import Fido, attrs as a

from astropy.coordinates import SkyCoord
from astropy import units as u

import warnings
warnings.filterwarnings("ignore")

result = Fido.search(a.Time('2014-09-19T05:59:00', '2014-09-19T06:01:00'), 
                     a.Instrument("aia"), a.Wavelength(171*u.angstrom), 
                     a.vso.Sample(12*u.second))
result