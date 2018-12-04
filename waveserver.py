#!/usr/bin/env python

STATIONS = [
    "NP.2704",
    "NP.2738",
    "NP.2784",
    "NP.8007",
    "NP.8011",
    "NP.8019",
    "NP.8021",
    "NP.8025",
    "NP.8026",
    "NP.8027",
    "NP.8028",
    "NP.8029",
    "NP.8030",
    "NP.8034",
    "NP.8036",
    "NP.8037",
    "NP.8038",
    "NP.8039",
    "NP.8041",
    "NP.8046",
    "NP.8047",
    "NP.8051",
    "NP.8052",
    "NP.ABBK",
    "NP.AHOU",
    "NP.ALUK",
    "NP.AMJG",
    "NP.ARTY",
    "NP.NIKO",
    "AK.CDVT",
    "AK.DAM1",
    "AK.DAM2",
    "AK.FA01",
    "AK.FA02",
    "AK.FA05",
    "AK.FA06",
    "AK.FA07",
    "AK.FA09",
    "AK.FA10",
    "AK.FA11",
    "AK.FA12",
    "AK.FIRE",
    "AK.GHO",
    "AK.GLI",
    "AK.HIN",
    "AK.PWL",
    "AK.TAPE",
    "AK.WAT1",
    "AK.K203",
    "AK.K204",
    "AK.K205",
    "AK.K208",
    "AK.K209",
    "AK.K210",
    "AK.K211",
    "AK.K212",
    "AK.K213",
    "AK.K214",
    "AK.K215",
    "AK.K216",
    "AK.K217",
    "AK.K218",
    "AK.K220",
    "AK.K221",
    "AK.K222",
    "AK.K223",
    "TA.N20K",
    ]
WAVESERVERS = [
    ("heli3.wr.usgs.gov", 16027), # Netquakes/Networms
    ("heli3.wr.usgs.gov", 16024), # AK network
]
from obspy import UTCDateTime, Stream
ot = UTCDateTime("2018-11-30T17:29:28")

from obspy.clients.earthworm import Client

stations_missing = list(STATIONS)
waveforms = Stream()

for server, port in WAVESERVERS:
    client = Client(server, port)
    response = client.get_availability(network="*", station="*", channel="*")
    
    stations_found = []
    for network, station, location, channel, start, end in response:
        label = ".".join([network, station])
        if label in STATIONS:
            stations_found.append(label)
            if label in stations_missing:
                stations_missing.remove(label)
            if UTCDateTime(start) < ot and UTCDateTime(end) > ot:
                st = client.get_waveforms(network=network, station=station, location=location, channel=channel, starttime=ot-30, endtime=ot+210)
                if st.count() > 0:
                    waveforms += st
                else:
                    print("No waveforms returned for station {label}.".format(label=label))
                    print("  starttime: {}, endtime: {}".format(start, end))
            else:
                print("No waveforms for station {label} on {server}:{port}.".format(
                    label=label, server=server, port=port))
                print("  starttime: {}, endtime: {}".format(start, end))
                                                                                        

    print("Waveserver: {}:{}".format(server, port))
    for st in stations_found:
        print("  {}".format(st))

print("Stations not found:")
for st in stations_missing:
    print("  {}".format(st))

if waveforms.count() > 0:
    filename = "data/ak20419010-Anchorage/waveserver.mseed"
    waveforms.write(filename, format="MSEED")
