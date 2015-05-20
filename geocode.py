"""
Author: Eric J. Ma

Purpose of Code: 
This code is a code snippet that shows how to search for GPS longitude 
and latitude coordinates based on a given set of search terms.

"""

from geopy.geocoders import GeoNames, GoogleV3

def get_geo_coordinates(search_terms):
	"""
	Input: 
	An iterable of terms to be searched. The search terms 
	should be country|locale, separated by a vertical pipe.
	Output: A list of tuples, where each tuple comprises the:
		- search term
		- returned location string based on the search term
		- longitude
		- latitude
	"""
	geocoder = GoogleV3(timeout=5)
	geodata = []
	for search_term in search_terms:
		new_search_term = search_term.replace('|', ' ')
		loc = geocoder.geocode(new_search_term)
		if loc is None:
			loc = geocoder.geocode(search_term.split('|')[0])
		lat = loc.latitude
		lon = loc.longitude
		data = (search_term, u'{0}'.format(loc), lon, lat)
		geodata.append(data)

	return geodata