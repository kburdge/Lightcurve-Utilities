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
