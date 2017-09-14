#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import os
import logging
import numpy

from obspy.clients.fdsn import Client

DEFAULTS = """
[waveformsapp]
event = nc72282711
name = South Napa

[fdsn.client]
data_centers = [NCEDC, IRIS]
debug = False

[fdsn.event]

[fdsn.station]
channel = HH?,HN?,EH?
maxradius = 0.5 ; degrees

[fdsn.dataselect]
location = *

[fdsn.event]

[stations]
#required_channels = [HHE,HHN,HHZ]
#priority_channels = [HN?,HH?]

[waveforms]
window_preevent = 10.0
window_total = 70.0

[map]
width_pixels = 4096
height_pixels = 4096
zoom_level = 11,
haxis_size = 0.05,
vaxis_size = 0.02,
time_window = 60.0,
show_axes = True,
show_labels = True,
mechanism_size = 100
mechanism_color = ltred

[record_section]
width = 7.5,
height = 10.0
time_window = 30.0
acc_per_km = 1.0
vel_per_km = 0.025
show_labels = True

[files]
event = event.xml
stations = stations_%%s.xml
waveforms_raw = waveforms_raw.mseed
waveforms_acc = waveforms_acc.p
waveforms_vel = waveforms_vel.p
maptiles = ~/scratch/images/tiles
"""

# ----------------------------------------------------------------------
def _config_get_list(list_string):
    """Convert list as string to list.

    :type list_string: list
    :param list_string: List as string.
    :returns: List of strings.
    """
    l = [f.strip() for f in list_string[1:-1].split(",")]
    return l

def _data_filename(params, pname, args=None, make_dir=False):
    """Construct relative path for file.
    
    :type params: ConfigParser
    :param params: Application parameters.
    :type pname: str
    :param pname: Name of parameter for file.
    :type args: tuple
    :param args: Tuple for arguments for substitution in name of parameter file.
    :type make_dir: bool
    :param make_dir: Create directory for file if True.
    """
    eventId = params.get("waveformsapp", "event")
    eventName = params.get("waveformsapp", "name").replace(" ","")
    eventDir = "%s-%s" % (eventId, eventName)

    if make_dir and not os.path.isdir(os.path.join("data", eventDir)):
        os.makedirs(os.path.join("data", eventDir))

    filename = params.get("files", pname) if args is None else params.get("files", pname) % args
    return os.path.join("data", eventDir, filename)
    

# ----------------------------------------------------------------------
class WaveformData(object):
    """Object for retrieving and processing waveform data.
    """

    def __init__(self, params, show_progress=True):
        """Constructor.

        :type params: ConfigParser
        :param params: Parameters for application.
        :type show_progress: bool
        :param show_progress: Show progress during execution.
        """
        self.params = params
        self.showProgress = show_progress

        self.event = None
        self.inventory = {}
        self.stream = None
        return

    def fetch_event(self):
        """Fetch event information from first data center in list of data centers.
        """
        import re
        
        dataCenter = _config_get_list(self.params.get("fdsn.client", "data_centers"))[0]
        
        if self.showProgress:
            print("Fetching event information from %s..." % dataCenter)

        client = Client(dataCenter, debug=self.params.getboolean("fdsn.client", "debug")) # Fetch from first data center
        eventid = self.params.get("waveformsapp", "event")
        kwargs = {
            "eventid": re.sub(r"\D", "", eventid),
        }
        kwargs.update(dict(self.params.items("fdsn.event")))
        catalog = client.get_events(**kwargs)
        assert(len(catalog.events) == 1)
        event = catalog.events[0]

        event.write(_data_filename(self.params, "event", make_dir=True), format="QUAKEML")

        self.event = event
        return

    def show_event(self):
        """Write earthquake information to stdout.
        """
        event = self._get_event()
        print(event)
        return
    
    def fetch_stations(self):
        """
        Get station information.
    
        @pre Must have fetched event.
        """
        import obspyutils.cisn
        
        if self.showProgress:
            print("Fetching station information from datacenter(s)...")

        hypocenter = self._select_hypocenter()

        self.inventory = {}
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            client = Client(dc, timeout=300, debug=self.params.getboolean("fdsn.client", "debug"))
            kwargs = {
                "startbefore": hypocenter.time,
                "endafter": hypocenter.time,
                "longitude": hypocenter.longitude,
                "latitude": hypocenter.latitude,
                "level": "response",
            }
            kwargs.update(dict(self.params.items("fdsn.station")))
            inventory = client.get_stations(**kwargs)
                
            obspyutils.cisn.remove_structures(inventory)

            # :TEMPORARY: remove PBO stations (no instrument response)
            for network in inventory.networks:
                if network.code == "PB":
                    network.stations = []

            inventory.write(_data_filename(self.params, "stations", dc), format="STATIONXML")

            if self.showProgress:
                numNetworks = len(inventory)
                numStations = 0
                for network in inventory.networks:
                    numStations += len(network.stations)
                print("    Retrieved %d total stations across %d networks from data center %s." % (numStations, numNetworks, dc))

            self.inventory[dc] = inventory
        return


    def show_stations(self):
        """Write stations information to stdout.
        """
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventory = self._get_inventory(dc)
            print(inventory)
        return
    
    def fetch_waveforms(self):
        if self.showProgress:
            print("Fetching recorded waveforms from data center(s)...")

        hypocenter = self._select_hypocenter()
        t1 = hypocenter.time - self.params.getfloat("waveforms", "window_preevent")
        t2 = t1 + self.params.getfloat("waveforms", "window_total")

        from math import ceil
        stream = None
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventory = self._get_inventory(dc)

            client = Client(dc, timeout=1200, debug=self.params.getboolean("fdsn.client", "debug"))
            bulk = []
            location = self.params.get("fdsn.dataselect", "location")
            for network in inventory.networks:
                for station in network.stations:
                    channels = self._select_channels(station)
                    if channels:
                        info = (network.code, station.code, location, channels, t1, t2)
                        bulk.append(info)
            maxSize = 200
            nbatches = int(ceil(len(bulk) / float(maxSize)))
            print "   Using %d batches to fetch waveforms for %d stations from data center %s..." % (nbatches, len(bulk), dc)
            for i in xrange(nbatches):
                istart = i*maxSize
                iend = min((i+1)*maxSize, len(bulk))
                streamB = client.get_waveforms_bulk(bulk[istart:iend])
                if stream is None:
                    stream = streamB
                else:
                    stream += streamB
            
        stream.write(_data_filename(self.params, "waveforms_raw"), format="MSEED")
        self.stream = stream
        return


    def process_waveforms(self):
        import obspyutils.rotate
        import obspyutils.metadata
        import pyproj
        
        if self.showProgress:
            print("Processing recorded waveforms...")

        hypocenter = self._select_hypocenter()

        inventory = None
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventoryDC = self._get_inventory(dc)
            if inventory is None:
                inventory = inventoryDC
            else:
                inventory += inventoryDC

        stream = self._get_stream()
        stream.attach_response(inventory)

        # Add metadata
        obspyutils.metadata.addLocation(inventory, stream)
        utmZone = int(math.floor(hypocenter.longitude + 180.0)/6.0) + 1
        proj = pyproj.Proj(proj="utm", zone=utmZone, ellps="WGS84")
        obspyutils.metadata.addAzimuthDist(stream, epicenter=(hypocenter.longitude, hypocenter.latitude), projection=proj)

        # Strong-motion acceleration traces
        smAcc = stream.select(channel="HN?")
        smAcc.remove_response(output="ACC")
        # Remove noise and perform baseline correction
        obspyutils.baselinecorrection.denoise(smAcc)
        # :TODO: add baseline correction
        smVel = smAcc.copy().integrate()

        # Broadband velocity traces
        bbVel = stream.select(channel="HH?")
        bbVel += stream.select(channel="EH?")
        bbVel.remove_response(output="VEL")
        bbAcc = bbVel.copy().differentiate()

        # Combine strong-motion and broadband streams
        self.acc = smAcc + bbAcc
        self.vel = smVel + bbVel
        
        # Rotate
        for s in [acc, vel]:
            obspyutils.rotate.toENZ(inventory, s)

        with open(_data_filename(self.params, "acc"), "w") as fout:
            cPickle.Pickler(fout, protocol=-1).dump(self.acc)
        with open(_data_filename(self.params, "vel"), "w") as fout:
            cPickle.Pickler(fout, protocol=-1).dump(self.vel)
        return

    def _get_event(self):
        """Read event if not already loaded.

        :returns: Event object
        """
        if self.event is None:
            import obspy.core.event
            self.event = obspy.core.event.read_events(_data_filename(self.params, "event"), format="QUAKEML")[0]
        return self.event

    def _get_inventory(self, dataCenter):
        """Read inventory information if not already loaded.

        :type dataCenter: str
        :param dataCenter: Data center associated with inventory.

        :returns: Inventory object
        """
        if not dataCenter in self.inventory:
            import obspy.core.inventory
            self.inventory[dataCenter] = obspy.core.inventory.read_inventory(_data_filename(self.params, "stations", dataCenter), format="STATIONXML")
        return self.inventory[dataCenter]
            
    def _get_stream(self):
        """Read stream if not already loaded.

        :returns: Waveform data stream
        """
        if self.stream is None:
            import obspy.core.stream
            self.stream = obspy.core.stream.read(self.params.get("files", "waveforms_raw"), format="MSEED")
        return self.stream

    def _select_channels(self, station):
        """Select preferred channels for station.

        :type station: obspy.core.inventory.station.Station
        :param station: Station information.
        """
        channels = station.channels
        if self.params.has_option("stations", "channel_required"):
            channels = self.params.get("stations", "channel_required")
        elif self.params.has_option("stations", "channel_priority"):
            channels = []
            for target in self.params.get("stations", "channel_priority").split(","):
                for channel in station.channels:
                    if channel.startswith(target):
                        channels.append(channel)
                if len(channels) > 0:
                    break
            channels = ",".join(channels)

        logging.getLogger(__name__).debug("Selected channels %s for station %s." % (channels, station.code,))
        return channels


    def _select_hypocenter(self):
        """Select hypocenter from preferred list of hypocenters.
        """
        event = self._get_event()
        hypocenter = event.preferred_origin()
        if hypocenter is None:
            hypocenter = event.origins[0]
        logging.getLogger(__name__).debug("Selected hypocenter %s." % hypocenter.method_id)
        return hypocenter

# ----------------------------------------------------------------------
class NoiseFigure(object):
    """
    """

# ----------------------------------------------------------------------
class WaveformsMap(object):
    """
    """

    def __init__(self, data):
        self.data = data
        return
    
    
    def plot(self):
        import basemap.Tiler
        import basemap.WaveformsMap

        paramsMap = self.params["map"]
        plotW = paramsMap["width_pixels"]
        plotH = paramsMap["height_pixels"]

        hypocenter = self._select_hypocenter()
        event = self._get_event()

        center = (hypocenter.longitude, hypocenter.latitude)
        domainH = 150.0 #2.5*self.params["stations"]["maxradius"]

        tiler = Tiler(center=center, width=plotW/float(plotH)*domainH, height=domainH, level=paramsMap["zoom_level"], base="ESRI_streetmap")
        tiler.initialize()
        tiler.tile("basemap.png")

        freqMin = self.params["waveforms"]["freq_min"]
        freqMax = self.params["waveforms"]["freq_max"]
        for field in ["acc","vel"]:
            stream = obspyutils.pickle.unpickle(self.params["waveforms"]["%s" % field])
            stream.rotate("NE->RT")
            if freqMax:
                if freqMin:
                    stream.filter("bandpass", freqmin=freqMin, freqmax=freqMax, corners=2, zerophase=True)
                else:
                    stream.filter("lowpass", freq=freqMax, corners=2, zerophase=True)

            elif freqMin:
                stream.filter("highpass", freq=freqMin, corners=2, zerophase=True)

            starttime = hypocenter.time
            endtime = starttime + paramsMap["time_window"]
            stream.trim(starttime, endtime, pad=True)

            for component in ["R","T","Z"]:
                streamC = stream.select(component=component)

                map = WaveformsMap(tiler, color="lightbg", fontsize=8)
                map.open(plotW, plotH, "basemap.png")

                from obspy.imaging.beachball import beach
                import obspyutils.event
                fm = obspyutils.event.first_motion(event)
                mechanism = (fm.nodal_planes.nodal_plane_1.strike.real,
                             fm.nodal_planes.nodal_plane_1.dip.real,
                             fm.nodal_planes.nodal_plane_1.rake.real)
                xy = tiler.geoToImage(numpy.array([[hypocenter.longitude, hypocenter.latitude]]))
                beachball = beach(mechanism, xy=xy[0], width=paramsMap["mechanism_size"], facecolor=paramsMap["mechanism_color"], axes=map.ax)
                map.ax.add_collection(beachball)
 
                map.plotStations(stream, showLabel=paramsMap["show_labels"])
                maxamp = streamC.max()
                if field == "acc":
                    maxlim = numpy.max(maxamp[maxamp < 2000.0])
                    ylabel = "Acceleration (m/s**2)"
                else:
                    maxlim = numpy.max(maxamp[maxamp < 2.0])
                    ylabel = "Velocity (m/s)"
                map.plotWaveforms(streamC, width=paramsMap["haxis_size"], height=paramsMap["vaxis_size"], showFirstAxes=True, ylim=[-maxlim,maxlim], label=ylabel)

                map.figure.savefig("waveformsmap_%s_%s.png" % (field, component))
    
# ----------------------------------------------------------------------
class RecordSection(object):
    """
    """

    def __init__(self, data):
        self.data = data
        return
    
    
    def plot(self):
        from basemap.Figure import Figure

        paramsRecSec = self.params["record_section"]

        hypocenter = self._select_hypocenter()

        for field in ["acc","vel"]:
            stream = obspyutils.pickle.unpickle(self.params["waveforms"]["%s" % field])
            stream.rotate("NE->RT")

            starttime = hypocenter.time
            endtime = starttime + paramsRecSec["time_window"]
            stream.trim(starttime, endtime, pad=True)

            for component in ["R","T","Z"]:
                streamC = stream.select(component=component)

                fig = Figure(color="lightbg", fontsize=10)
                fig.open(width=paramsRecSec["width"], height=paramsRecSec["height"], margins=((0.5, 0, 0.1), (0.5, 0, 0.1)))
                ax = fig.axes(1, 1, 1, 1)
                
                vscale = paramsRecSec["%s_per_km" % field]

                for trace in streamC.traces:
                    distKm = trace.stats.distance / 1.0e+3
                    if numpy.max(numpy.abs(trace.data/vscale)) < 5.0:
                        ax.plot(trace.times(), distKm+trace.data/vscale, linewidth=0.5, color="blue")

                        if paramsRecSec["show_labels"]:
                            label = "%s.%s" % (trace.stats.network, trace.stats.station)
                            ax.text(numpy.min(trace.times()), distKm, label, horizontalalignment="left", verticalalignment="bottom", fontsize=6, color="orange")

                ax.set_xlabel("Time (s)")
                ax.set_ylabel("Epicentral Distance (km)")
                ax.set_xlim(0.0, paramsRecSec["time_window"])
                fig.figure.savefig("recordsection_%s_%s.pdf" % (field, component))
    

# ----------------------------------------------------------------------
class WaveformsApp(object):
    """
    Plot record section and map with acceleration and velocity waveforms.
    """
    
    def __init__(self, show_progress=True):
        """Constructor.

        :type show_progress: bool
        :param show_progress: Show progress during execution.
        """
        self.showProgress = show_progress
        return

    def initialize(self, config_filenames):
        """Set parameters from config file and DEFAULTS.

        :type config_filename: str
        :param config_filename: Name of configuration (INI) file with parameters.
        """
        import ConfigParser
        import io
        config = ConfigParser.SafeConfigParser()
        config.readfp(io.BytesIO(DEFAULTS))
        for filename in config_filenames.split(","):
            if self.showProgress:
                print("Fetching parameters from %s..." % filename)
            config.read(filename)

        self.params = config
        return
    
    def show_parameters(self):
        """Write parameters to stdout.
        """
        import sys
        self.params.write(sys.stdout)
        return

    def event_name(self):
        """Get event name from event id and earthquake name/location.

        :return: Event name
        """
        eventid = self.params.get("waveformsapp", "event")
        eventname = self.params.get("waveformsapp", "name").replace(" ","").lower()
        return "%s-%s" % (eventid, eventname)
        

# ======================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", action="store", dest="config", required=True)
    parser.add_argument("--show-parameters", action="store_true", dest="show_parameters")
    parser.add_argument("--fetch-event", action="store_true", dest="fetch_event")
    parser.add_argument("--show-event", action="store_true", dest="show_event")
    parser.add_argument("--fetch-stations", action="store_true", dest="fetch_stations")
    parser.add_argument("--show-stations", action="store_true", dest="show_stations")
    parser.add_argument("--fetch-waveforms", action="store_true", dest="fetch_waveforms")
    parser.add_argument("--process-waveforms", action="store_true", dest="process_waveforms")
    parser.add_argument("--plot-section", action="store_true", dest="plot_section")
    parser.add_argument("--plot-noise", action="store_true", dest="plot_noise")
    parser.add_argument("--plot-map", action="store_true", dest="plot_map")
    parser.add_argument("--all", action="store_true", dest="all")
    parser.add_argument("--quiet", action="store_false", dest="show_progress", default=True)
    parser.add_argument("--debug", action="store_true", dest="debug")
    args = parser.parse_args()

    app = WaveformsApp(args.show_progress)
    app.initialize(args.config)

    logLevel = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=logLevel, filename=app.event_name()+".log")

    data = WaveformData(app.params, args.show_progress)
    
    if args.show_parameters or args.all:
        app.show_parameters()
    
    if args.fetch_event or args.all:
        data.fetch_event()

    if args.show_event or args.all:
        data.show_event()

    if args.fetch_stations or args.all:
        data.fetch_stations()

    if args.show_stations or args.all:
        data.show_stations()

    if args.fetch_waveforms or args.all:
        data.fetch_waveforms()

    if args.process_waveforms or args.all:
        data.process_waveforms()

    if args.plot_noise or args.all:
        plotter = NoiseFigure(data)
        plotter.plot()
        
    if args.plot_map or args.all:
        plotter = WaveformsMap(data)
        plotter.plot()
    
    if args.plot_section or args.all:
        plotter = RecordSection(data)
        plotter.plot()
    

# End of file
