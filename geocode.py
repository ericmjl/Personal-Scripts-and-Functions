"""
Author: Eric J. Ma

Purpose of Code: 
This is a function that gets the longitude and latitude of a search term. 
"""

from geopy.geocoders import GeoNames, GoogleV3

def get_lon_lat_coordinates(search_term, engine):
	"""
	Gets the longitude and latitude of a specified search term.
	
	Parameters:
	===========

	-	search_term:	string. The name of the place for which you want to search.
	-	engine:			string. Currently supports:
						- 'geonames'
						- 'google'

	Returns:
	========
	-	(lon, lat)		tuple of floats. The longitude and latitude coordinates.

	"""
	gc = {'geonames':GeoNames(country_bias='USA', username='ericmjl'),
		  'google':GoogleV3(timeout=5)}
	geocoder = gc[engine]
	loc = geocoder.geocode(search_term)
	lat = loc.latitude
	lon = loc.longitude

	return lon, lat