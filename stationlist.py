#!/usr/bin/env python

import numpy


FILENAME = "data/ak20419010-Anchorage/stationlist_ak.json"

EPICENTER = (-149.923, 61.323)

EARTH_MEAN_RADIUS_M = 6371.0e+3
DEG_TO_RAD = numpy.pi / 180.0


def distance(refLon, refLat, ptsLon, ptsLat):
    """Get great circle distance in meters from reference point to points.

    Source: https://en.wikipedia.org/wiki/Great-circle_distance

    :type refLon: float
    :param refLon: Longitude of reference point in degrees.

    :type refLat: float
    :param refLat: Latitude of reference point in degrees.

    :type ptsLon: Numpy array
    :param ptsLon: Longitude of points in degrees.

    :type ptsLat: Numpy array
    :param ptsLat: Latitude of points in degrees.
    """
    refLonR = refLon * DEG_TO_RAD
    refLatR = refLat * DEG_TO_RAD
    ptsLonR = ptsLon * DEG_TO_RAD
    ptsLatR = ptsLat * DEG_TO_RAD

    p = numpy.sin(0.5*(ptsLatR-refLatR))**2 \
        + numpy.cos(refLatR)*numpy.cos(ptsLatR)*numpy.sin(0.5*(ptsLonR-refLonR))**2
    return EARTH_MEAN_RADIUS_M * 2.0*numpy.arcsin(p**0.5)

stations = [
    ("NP.2716","Hilton, Anchorage",-149.8921, 61.2193),
    ("NP.2750","Alaska Regional Hospital, Anchorage",-149.8282,61.2107),
    ("NP.8016","BP Building, Anchorage",-149.8645,61.1922),
    ("NP.8040","Atwood Building, Anchorage",-149.8933,61.2135),
    ("NP.8042","Frontier Building, Anchorage",-149.8837,61.1876),
    ("NP.8043","Port Access Bridge, Anchorage",-149.8847,61.2218),
    ("NP.8045","VAMC, Anchorage",-149.7439,61.2332),
    ("NP.2704","Old Federal Bld",-149.8941,61.2191),
    ("NP.2738","ADOT Maint Sta, Cantwell, AK",-148.8849,63.3888),
    ("NP.2784","Valdez City Hall, Valdez, AK",-146.3547,61.1303),
    ("NP.8019","Alaska Geol Mat, Eagle River, AK",-149.5411,61.3493),
    ("NP.8034","Valdez Civit Ctr, Valdez, AK",-146.3567,61.1263),
    ("NP.8039","FS 07, Anchorage, AK",-149.9512,61.1416),
    ("NP.8046","VAMC free field, Anchorage, AK",-149.7435,61.2338),
    ("NP.8047","USGS ESC, Anchorage, AK",-149.8020,61.1885),
    ("NP.NIKO","Nikolaevsk School, Nikolaevsk, AK",-151.6127,59.8113),
    ("AK.FIRE","Fire Island, AK",-150.2164,61.1426),
    ("AK.GHO","Gloryhole, AK",-148.9260,61.7710),
    ("AK.HIN","Hinchinbrook, AK",-146.5035,60.3960),
    ("AK.PWL","Port Wells, AK",-148.3334,60.8584),
    ("TA.N20K","Mount Spurr, AK",-152.2089,61.2001),
    ("NP.8007","Anchorage Intl Airpport, Anchorage, AK",-149.9967,61.1824),
]
for code,name,lon,lat in stations:
    dist = 1.0e-3*distance(EPICENTER[0], EPICENTER[1], lon, lat)
    print("{code} {dist:.2f}".format(code=code, dist=dist))

import json
with open(FILENAME, "r") as fin:
    data = json.load(fin)

    for feature in data["features"]:
        props = feature["properties"]
    
        network = props["network"]
        station_code = props["code"]
        station_name = props["name"]
        dist_km = props["distance"]
        channels = ",".join([channel["name"] for channel in props["channels"]])

        coords = feature["geometry"]["coordinates"]
        longitude = coords[0]
        latitude = coords[1]

        epicentral_dist = 1.0e-3 * distance(EPICENTER[0], EPICENTER[1], longitude, latitude)
        
        if channels.startswith("HN"):
            print("{network}.{station_code}|{name}|{longitude:.4f}|{latitude:.4f}|{channels}|{ep_dist:.2f}|{sm_dist:.2f}".format(
                network=network, station_code=station_code, name=station_name, longitude=longitude, latitude=latitude,
                channels=channels, ep_dist=epicentral_dist, sm_dist=dist_km))
