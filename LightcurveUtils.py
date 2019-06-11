#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 09:16:09 2019

@author: kburdge
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5, BarycentricTrueEcliptic  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u


"Function for converting a series of timestamps to Barycentric Julian Date format in Barycentric Dynamical time"

def BJDConvert(times, RA, Dec, date_format='mjd', telescope='Palomar'):

	t = Time(times,format=date_format,scale='utc')

	t2=t.tdb

	c = SkyCoord(RA,Dec, unit="deg")

	d=c.transform_to(BarycentricTrueEcliptic)

	Observatory=EarthLocation.of_site(telescope)

	delta=t2.light_travel_time(c,kind='barycentric',location=Observatory)

	BJD_TDB=t2+delta

	return BJD_TDB

"Function which returns phases corresponding to timestamps in a lightcurve given a period P, period derivative Pdot, and reference epoch t0"

def PhaseFold(times, P, Pdot=0, t0=0):


	"If no reference epoch is supplied, reference epoch is set to earliest time in lightcurve"
	if t0==0:
		times=times-np.min(times)
	else:
		times=times-t0
	phases=((times-1/2*Pdot/P*(times)**2) % P)/P

	return phases
