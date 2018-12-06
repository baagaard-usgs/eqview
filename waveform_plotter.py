#!/usr/bin/env python
# ======================================================================
# Brad T. Aagaard
# U.S. Geological Survey
# ======================================================================

import os
import logging
import numpy

from obspy.clients.fdsn import Client

DEFAULTS = """
[waveformsapp]
event = nc72282711
name = South Napa

[fdsn.client]
data_centers = [USGS, NCEDC, IRIS]
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

[processing.noise]
remove_bg = True
zero_coarse_levels = 1
zero_fine_levels = 1
preevent_threshold_reduction = 2.0
store_noise = False
store_orig = False
wavelet = coif4

[processing.filter]
zero_phase = True
num_corners = 2
freq_min = 0.2
freq_max = 25.0

[maps]
tiler = cartopy_extra_tiles.esri_tiles.ESRI
tiler_style = streetmap
tiler_cache_dir = ~/data_scratch/images/tiles
width_in = 15.0
height_in = 15.0
zoom_level = 8

[waveforms_map]
haxis_size = 0.05
vaxis_size = 0.02
show_axes = True
show_labels = False
mechanism_size = 10
station_size = 6

[noise_figure]
height = 7.5
width = 20.0
margins = ((0.6, 0.5, 0.2), (0.5, 0.6, 0.7))


[record_section]
width_in = 15.0
height_in = 15.0
margins = ((0.8, 0.5, 0.2), (0.5, 0.6, 0.7))

time_window = 60.0
acc_per_km = 0.3
vel_per_km = 0.025
show_labels = True
show_wavespeeds = True
vp_kmps = 5.6
vs_kmps = 3.4

[files]
event = event.xml
stations = stations_%%s.xml
waveforms_raw = waveforms_raw.mseed
waveforms_acc = waveforms_acc.p
waveforms_vel = waveforms_vel.p
plots = plots
maptiles = ~/data_scratch/images/tiles
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

        self.accSM = None
        self.velBB = None
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
            "eventid": eventid if dataCenter == "USGS" else re.sub(r"\D", "", eventid),
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
            if dc == "USGS":
                continue
            
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
                print("    Retrieved {} total stations across {} networks from data center {}.".format(numStations, numNetworks, dc))

            self.inventory[dc] = inventory
        return


    def show_stations(self):
        """Write stations information to stdout.
        """
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventory = self._get_inventory(dc)
            if inventory is None:
                continue

            for network in inventory.networks:
                for station in network:
                    channels = ",".join(["{}.{}".format(channel.location_code, channel.code) for channel in station.channels])
                    print("{network.code}.{station.code} {station.site.name} {station.longitude:.4f} {station.latitude:.4f} {channels}".format(
                        network=network, station=station, channels=channels))
        return
    
    def fetch_waveforms(self):
        if self.showProgress:
            print("Fetching recorded waveforms from data center(s)...")

        hypocenter = self._select_hypocenter()
        t1 = hypocenter.time - self.params.getfloat("waveforms", "window_preevent")
        t2 = t1 + self.params.getfloat("waveforms", "window_total")

        from math import ceil
        from obspy import Stream
        
        stream = Stream()
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventory = self._get_inventory(dc)
            if inventory is None:
                continue

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
            print("   Using {} batches to fetch waveforms for {} stations from data center {}...".format(nbatches, len(bulk), dc))
            for i in xrange(nbatches):
                istart = i*maxSize
                iend = min((i+1)*maxSize, len(bulk))
                stream += client.get_waveforms_bulk(bulk[istart:iend])
            
        stream.write(_data_filename(self.params, "waveforms_raw"), format="MSEED")
        self.stream = stream
        return


    def process_waveforms(self):
        import obspyutils.rotate
        import obspyutils.metadata
        import obspyutils.noise
        import obspyutils.baseline
        import pyproj
        import cPickle
        import math
        
        if self.showProgress:
            print("Processing recorded waveforms...")

        hypocenter = self._select_hypocenter()

        inventory = None
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventoryDC = self._get_inventory(dc)
            if inventoryDC is None:
                continue
            
            if inventory is None:
                inventory = inventoryDC
            else:
                inventory += inventoryDC

        stream = self._get_raw_stream()
        for bad_station in _config_get_list(self.params.get("stations", "blacklist")):
            net, st, loc, ch = bad_station.split(".")
            for tr in stream.select(network=net, station=st, location=loc, channel=ch):
                stream.remove(tr)
                                       
        stream.attach_response(inventory)

        # Add metadata
        obspyutils.metadata.addLocation(inventory, stream)
        utmZone = int(math.floor(hypocenter.longitude + 180.0)/6.0) + 1
        proj = pyproj.Proj(proj="utm", zone=utmZone, ellps="WGS84")
        obspyutils.metadata.addAzimuthDist(stream, epicenter=(hypocenter.longitude, hypocenter.latitude), projection=proj)

        # Strong-motion acceleration traces
        sSM = stream.select(channel="HN?")
        sSM.remove_response(output="ACC")
        # Remove noise and perform baseline correction
        wavelet = self.params.get("processing.noise", "wavelet")
        removeBg = self.params.getboolean("processing.noise", "remove_bg")
        preWindow = self.params.getfloat("waveforms", "window_preevent")
        preThreshold = self.params.getfloat("processing.noise", "preevent_threshold_reduction")
        zeroCoarse = self.params.getint("processing.noise", "zero_coarse_levels")
        zeroFine = self.params.getint("processing.noise", "zero_fine_levels")
        storeNoise = self.params.getboolean("processing.noise", "store_noise")
        storeOrig = self.params.getboolean("processing.noise", "store_orig")
        accSM = obspyutils.noise.denoise(sSM, wavelet, removeBg, zeroCoarse, zeroFine, preWindow, preThreshold, storeOrig, storeNoise)

        freqMin = self.params.getfloat("processing.filter", "freq_min")
        freqMax = self.params.getfloat("processing.filter", "freq_max")
        numCorners = self.params.getint("processing.filter", "num_corners")
        zeroPhase = self.params.getboolean("processing.filter", "zero_phase")
        if "orig" in accSM:
            accSM["bandpass"] = accSM["orig"].copy()
            accSM["bandpass"].filter("bandpass", freqmin=freqMin, freqmax=freqMax, corners=numCorners, zerophase=zeroPhase)
        
        obspyutils.baseline.correction_constant(accSM["data"])
        if "orig" in accSM:
            obspyutils.baseline.correction_constant(accSM["orig"])
        if "bandpass" in accSM:
            obspyutils.baseline.correction_constant(accSM["bandpass"])
        
        # Broadband velocity traces
        sVel = stream.select(channel="HH?")
        sVel += stream.select(channel="EH?")
        sVel.remove_response(output="VEL")
        velBB = obspyutils.noise.denoise(sVel, wavelet, removeBg, zeroCoarse, zeroFine, preWindow, preThreshold, storeOrig, storeNoise)
        
        # Rotate
        for s in accSM.values() + velBB.values():
            s = obspyutils.rotate.toENZ(inventory, s)

        self.accSM = accSM
        with open(_data_filename(self.params, "waveforms_acc"), "w") as fout:
            cPickle.Pickler(fout, protocol=-1).dump(accSM)

        self.velBB = velBB
        with open(_data_filename(self.params, "waveforms_vel"), "w") as fout:
            cPickle.Pickler(fout, protocol=-1).dump(velBB)
            
        return

    def load_processed_waveforms(self):
        """Load processed data.
        """
        import cPickle
        if self.showProgress:
            print("Loading processed waveforms...")            

        if self.accSM is None:
            with open(_data_filename(self.params, "waveforms_acc"), "r") as fin:
                self.accSM = cPickle.Unpickler(fin).load()
        if self.velBB is None:
            with open(_data_filename(self.params, "waveforms_vel"), "r") as fin:
                self.velBB = cPickle.Unpickler(fin).load()

        return
    
    def load_event(self):
        """Load event.
        """
        if self.showProgress:
            print("Loading event...")            

        self.event = self._get_event()
        self.hypocenter = self._select_hypocenter()
        return

    def load_stations(self):
        """Load stations. Use raw waveforms to get station locations if
        available, otherwise use list of stations.
        """
        import obspyutils.metadata
        
        if self.showProgress:
            print("Loading stations...")            

        inventory = None
        for dc in _config_get_list(self.params.get("fdsn.client", "data_centers")):
            inventoryDC = self._get_inventory(dc)
            if inventoryDC is None:
                continue
            
            if inventory is None:
                inventory = inventoryDC
            else:
                inventory += inventoryDC
            
        if os.path.isfile(_data_filename(self.params, "waveforms_raw")):
            stream = self._get_raw_stream()
            self.stationsSM = stream.select(channel="HN?")
            self.stationsBB = stream.select(channel="HH?")
            obspyutils.metadata.addLocation(inventory, self.stationsSM)
            obspyutils.metadata.addLocation(inventory, self.stationsBB)
        else:
            self.stationsSM = inventory.select(channel="HN?")
            self.stationsBB = inventory.select(channel="HH?")
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
        if dataCenter == "USGS":
            return
        
        if not dataCenter in self.inventory:
            import obspy.core.inventory
            self.inventory[dataCenter] = obspy.core.inventory.read_inventory(_data_filename(self.params, "stations", dataCenter), format="STATIONXML")
        return self.inventory[dataCenter]
            
    def _get_raw_stream(self):
        """Read raw stream if not already loaded.

        :returns: Waveform data stream
        """
        if self.stream is None:
            import obspy.core.stream
            self.stream = obspy.core.stream.read(_data_filename(self.params, "waveforms_raw"), format="MSEED")
        return self.stream

    def _select_channels(self, station):
        """Select preferred channels for station.

        :type station: obspy.core.inventory.station.Station
        :param station: Station information.
        """
        channels = ",".join([channel.code for channel in station.channels])
        if self.params.has_option("stations", "channel_required"):
            channels = self.params.get("stations", "channel_required")
            for target in channels.split(","):
                found = False
                for channel in station.channels:
                    if channel.code == target:
                        found = True
                        break
                if not found:
                    channels = None
                    logging.getLogger(__name__).debug("Missing required channel %s for station %s." % (target, station.code,))
                    break
        elif self.params.has_option("stations", "channel_priority"):
            channels = []
            for target in self.params.get("stations", "channel_priority").split(","):
                for channel in station.channels:
                    if channel.code.startswith(target):
                        channels.append(channel.code)
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
    Rows are original, denoised, and noise.
    Columns are coefficients, acceleration, velocity, and displacement.
    """

    ROWS = ["Original", "Signal", "Noise"]
    COLS = ["coefs", "acc", "vel", "disp"]
    
    def __init__(self, params, showProgress):
        """
        """
        self.params = params
        self.showProgress = showProgress
        return


    def plot(self, data):
        """
        """
        from ast import literal_eval
        import pywt
        import sys

        import matplotlib.pyplot as pyplot

        import obspyutils.baseline
        import matplotlib_extras.colors
        import matplotlib_extras
        
        if self.showProgress:
            sys.stdout.write("Plotting noise figures...")
        
        w = self.params.getfloat("noise_figure", "width")
        h = self.params.getfloat("noise_figure", "height")
        self.figure = pyplot.figure(figsize=(w,h))

        self._setupSubplots()
        self.figure.figure.canvas.draw()

        originTime = data._select_hypocenter().time
        preWindow = self.params.getfloat("waveforms", "window_preevent")

        numTraces = len(data.accSM["data"]) + len(data.velBB["data"])
        iTrace = 0
        for trAcc,trAccOrig,trAccNoise in zip(data.accSM["data"],data.accSM["orig"],data.accSM["noise"]):
            for trA,trB in ((trAcc,trAccOrig),(trAcc,trAccNoise),):
                assert(trA.stats.network == trB.stats.network)
                assert(trA.stats.station == trB.stats.station)
                assert(trA.stats.channel == trB.stats.channel)
            
            info = "%s.%s.%s\nS/N %.1f" % (trAcc.stats.network, trAcc.stats.station, trAcc.stats.channel, trAcc.StoN,)
            self.figure.figure.suptitle(info, fontweight='bold')

            dataRows = {
                "Original": trAccOrig,
                "Signal": trAcc,
                "Noise": trAccNoise,
            }

            for row in NoiseFigure.ROWS:
                # Coefficients
                coefs = pywt.wavedec(dataRows[row], self.params.get("processing.noise", "wavelet"), mode="zero")
                cArray,cSlices = pywt.coeffs_to_array(coefs)
                self._updatePlot(row+"_coefs", None, cArray)

                acc = dataRows[row]

                t = acc.times(reftime=originTime)
                vel, disp = obspyutils.baseline.integrate_acc(acc)
                self._updatePlot(row+"_acc", t, acc)
                self._updatePlot(row+"_vel", t, vel)
                self._updatePlot(row+"_disp", t, disp)

            plotsDir = os.path.join(_data_filename(self.params, "plots"))
            if not os.path.isdir(plotsDir):
                os.makedirs(plotsDir)
            filename = "%s.%s.%s.png" % (trAcc.stats.network, trAcc.stats.station, trAcc.stats.channel,)
            self.figure.figure.savefig(os.path.join(plotsDir, filename))

            if self.showProgress:
                sys.stdout.write("\rPlotting noise figures...%d%%" % (((iTrace+1)*100)/numTraces))
                sys.stdout.flush()
            iTrace += 1

        for trVel,trVelOrig,trVelNoise in zip(data.velBB["data"],data.velBB["orig"],data.velBB["noise"]):
            for trA,trB in ((trVel,trVelOrig),(trVel,trVelNoise),):
                assert(trA.stats.network == trB.stats.network)
                assert(trA.stats.station == trB.stats.station)
                assert(trA.stats.channel == trB.stats.channel)
            
            info = "%s.%s.%s\nS/N %.1f" % (trVel.stats.network, trVel.stats.station, trVel.stats.channel, trVel.StoN,)
            self.figure.figure.suptitle(info, fontweight='bold')

            dataRows = {
                "Original": trVelOrig,
                "Signal": trVel,
                "Noise": trVelNoise,
            }

            for row in NoiseFigure.ROWS:
                # Coefficients
                coefs = pywt.wavedec(dataRows[row], self.params.get("processing.noise", "wavelet"), mode="zero")
                cArray,cSlices = pywt.coeffs_to_array(coefs)
                self._updatePlot(row+"_coefs", None, cArray)

                vel = dataRows[row]

                t = vel.times(reftime=originTime)
                disp,dispI = obspyutils.baseline.integrate_acc(vel)
                self._updatePlot(row+"_acc", [], [])
                self._updatePlot(row+"_vel", t, vel)
                self._updatePlot(row+"_disp", t, disp)

            plotsDir = os.path.join(_data_filename(self.params, "plots"))
            if not os.path.isdir(plotsDir):
                os.makedirs(plotsDir)
            filename = "%s.%s.%s.png" % (trVel.stats.network, trVel.stats.station, trVel.stats.channel,)
            self.figure.figure.savefig(os.path.join(plotsDir, filename))

            if self.showProgress:
                sys.stdout.write("\rPlotting noise figures...%d%%" % (((iTrace+1)*100)/numTraces))
                sys.stdout.flush()
            iTrace += 1
            
        if self.showProgress:
            sys.stdout.write("\n")
            
        return

    def _setupSubplots(self):
        """
        """
        TITLES = {
            "coefs": "Coefficients (all levels)",
            "acc": "Acceleration (m/s**2)",
            "vel": "Velocity (m/s)",
            "disp": "Displacement (m)",
        }
        XLABELS = {
            "coefs": "",
            "acc": "Time (s)",
            "vel": "Time (s)",
            "disp": "Time (s)",
        }
        nrows = len(NoiseFigure.ROWS)
        ncols = len(NoiseFigure.COLS)

        margins = literal_eval(self.params.get("noise_figure", "margins"))
        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=margins)
        
        self.axes = {}
        for irow,row in enumerate(NoiseFigure.ROWS):
            for icol,col in enumerate(NoiseFigure.COLS):
                ax = pyplot.axes(rectFactory.rect(), nrows=nrows, ncols=ncols, row=irow+1, col=icol+1)
                line, = ax.plot([], [], 'r-', lw=0.5)
                ax.autoscale(enable=True, axis="both", tight=True)
                if irow == 0:
                    ax.set_title(TITLES[col])
                if irow == nrows-1:
                    ax.set_xlabel(XLABELS[col])
                if icol == 0:
                    pos = ax.get_position()
                    ax.text(pos.xmin, pos.ymax+0.02, row, fontweight='bold', transform=self.figure.figure.transFigure, ha="right")
                    ax.set_xticks([])
                self.axes[row+"_"+col] = (ax,line)
        return

    def _updatePlot(self, axiskey, x, y):
        ax, line = self.axes[axiskey]
        if x is not None:
            line.set_xdata(x)
        else:
            line.set_xdata(numpy.linspace(start=0,stop=y.shape[-1],num=y.shape[-1]))
        line.set_ydata(y)
        ax.relim()
        ax.autoscale_view()
        ax.draw_artist(line)
        return
    
# ----------------------------------------------------------------------
class WaveformsMap(object):
    """
    """

    def __init__(self, params, showProgress):
        self.params = params
        self.showProgress = showProgress
        return
    
    
    def plot(self, data):
        import sys

        if self.showProgress:
            sys.stdout.write("Plotting station map...")

        hypocenter = data.hypocenter
            
        plotW = self.params.getfloat("map", "width_pixels")
        plotH = self.params.getfloat("map", "height_pixels")
        center = (hypocenter.longitude, hypocenter.latitude)
        domainH = 2.5*110.0*self.params.getfloat("fdsn.station", "maxradius")
        zoomLevel = self.params.getint("map", "zoom_level")
        baseImage = self.params.get("map", "base_image")
        
        tiler = Tiler(center=center, width=plotW/float(plotH)*domainH, height=domainH, level=zoomLevel, base=baseImage)
        tiler.initialize()
        tiler.tile("basemap.png", cacheDir=os.path.expanduser(self.params.get("files", "maptiles")))

        plotsDir = os.path.expanduser(os.path.join(_data_filename(self.params, "plots")))
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)

        # Integrate accSM to velocity
        # Differentiate velBB to acceleration
    
        for field in ["acc","vel"]:
            starttime = hypocenter.time
            endtime = starttime + paramsMap["time_window"] # FIX THIS
            stream.trim(starttime, endtime, pad=True)

            for component in ["R","T","Z"]:
                streamC = stream.select(component=component) # FIX THIS

                map = WaveformsMap(tiler, color="lightbg", fontsize=8)
                map.open(plotW, plotH, "basemap.png")

                #from obspy.imaging.beachball import beach
                #import obspyutils.event
                #fm = obspyutils.event.first_motion(event)
                #mechanism = (fm.nodal_planes.nodal_plane_1.strike.real,
                #             fm.nodal_planes.nodal_plane_1.dip.real,
                #             fm.nodal_planes.nodal_plane_1.rake.real)
                #xy = tiler.geoToImage(numpy.array([[hypocenter.longitude, hypocenter.latitude]]))
                #beachball = beach(mechanism, xy=xy[0], width=paramsMap["mechanism_size"], facecolor=paramsMap["mechanism_color"], axes=map.ax)
                #map.ax.add_collection(beachball)
                map.plotEpicenter((hypocenter.longitude, hypocenter.latitude), size=self.params.getfloat("map", "mechanism_size"))
                map.plotStations(data.stationsSM, color="c_blue", marker="^", showLabels=True)
                map.plotStations(data.stationsBB, color="c_ltgreen", marker="v", showLabels=True)

                maxamp = streamC.max()
                if field == "acc":
                    maxlim = numpy.max(maxamp[maxamp < 2000.0])
                    ylabel = "Acceleration (m/s**2)"
                else:
                    maxlim = numpy.max(maxamp[maxamp < 2.0])
                    ylabel = "Velocity (m/s)"
                map.plotWaveforms(streamC, width=haxisSize, height=vaxisSize, showFirstAxes=True, ylim=[-maxlim,maxlim], label=ylabel)

                map.figure.savefig(os.path.join(plotsDir, "waveformsmap_%s_%s.png" % (field, component,)))
    
# ----------------------------------------------------------------------
class StationsMap(object):
    """
    """
    from cartopy import crs
    wgs84CRS = crs.Geodetic()
    
    def __init__(self, config, showProgress):
        self.config = config
        self.showProgress = showProgress
        return
    
    
    def plot(self, data):
        import sys
        from importlib import import_module

        import matplotlib.pyplot as pyplot

        from cartopy_extra_tiles import cached_tiler
        import matplotlib_extras
        
        if self.showProgress:
            sys.stdout.write("Plotting station map...")

        hypocenter = data.hypocenter


        tilerPath = self.config.get("maps", "tiler").split(".")
        tilerObj = getattr(import_module(".".join(tilerPath[:-1])), tilerPath[-1])
        tilerStyle = self.config.get("maps", "tiler_style")
        tilerZoom = self.config.getint("maps", "zoom_level")
        tilesDir = self.config.get("maps", "tiler_cache_dir")
        tiler = cached_tiler.CachedTiler(tilerObj(desired_tile_form="L", style=tilerStyle), cache_dir=tilesDir)

        figWidthIn = self.config.getfloat("maps", "width_in")
        figHeightIn = self.config.getfloat("maps", "height_in")
        figure = pyplot.figure(figsize=(figWidthIn, figHeightIn))
        matplotlib_extras.colors.add_general()

        rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=((0, 0, 0), (0, 0, 0)))
        ax = pyplot.axes(rectFactory.rect(), projection=tiler.crs)
        lon0 = hypocenter.longitude
        lat0 = hypocenter.latitude
        ax.set_extent([lon0-8.0, lon0+8.0, lat0-4.0, lat0+4.0])
        ax.add_image(tiler, tilerZoom, zorder=0, cmap="gray", interpolation="bessel")

        # Epicenter
        msize = self.config.getfloat("waveforms_map", "mechanism_size")
        ax.scatter(hypocenter.longitude, hypocenter.latitude, transform=self.wgs84CRS, marker="*", c="c_ltred", edgecolors="black", s=msize**2, zorder=4)

        # Stations
        self._plot_stations(ax, data.stationsSM, color="c_blue", marker="^")
        self._plot_stations(ax, data.stationsBB, color="c_ltgreen", marker="v")

        plotsDir = os.path.expanduser(os.path.join(_data_filename(self.config, "plots")))
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)
        figure.savefig(os.path.join(plotsDir, "stations_map.jpg"))
        return

    def _plot_stations(self, ax, stations, color, marker, show_labels=True, zorder=1):
        if 0 == stations.count():
            return
        
        lonlat = {}
        for tr in stations.traces:
            key = "%s.%s" % (tr.stats.network, tr.stats.station,)
            if not key in lonlat:
                lonlat[key] = (tr.stats.longitude, tr.stats.latitude,)
        coords = numpy.array(lonlat.values())
        msize = self.config.getfloat("waveforms_map", "station_size")
        ax.scatter(coords[:,0], coords[:,1], transform=self.wgs84CRS, marker="^", c=color, edgecolors="black", s=msize**2, alpha=0.7, zorder=zorder)

        if show_labels:
            for label,ll in zip(lonlat.keys(), lonlat.values()):
                ax.text(ll[0], ll[1], label+" ", transform=self.wgs84CRS, ha='right', va='center', color="black")
        
        return
    
# ----------------------------------------------------------------------
class RecordSection(object):
    """
    """

    def __init__(self, config, showProgress):
        self.config = config
        self.showProgress = showProgress
        return
    
    
    def plot(self, data):
        from ast import literal_eval

        import matplotlib.pyplot as pyplot

        import matplotlib_extras.colors
        import matplotlib_extras

        hypocenter = data.hypocenter

        # Acceleration
        field = "acc"
        stream = data.accSM["data"]
        tlength = self.config.getfloat("record_section", "time_window")
        starttime = hypocenter.time
        endtime = starttime + tlength
        stream.trim(starttime, endtime, pad=True, fill_value=0.0)
        stream.rotate("NE->RT")

        plotsDir = os.path.expanduser(os.path.join(_data_filename(self.config, "plots")))
        if not os.path.isdir(plotsDir):
            os.makedirs(plotsDir)

        for component in ["R","T","Z"]:
            streamC = stream.select(component=component)

            plotW = self.config.getfloat("record_section", "width_in")
            plotH = self.config.getfloat("record_section", "height_in")
            figure = pyplot.figure(figsize=(plotW,plotH))
            matplotlib_extras.colors.add_general()

            margins = literal_eval(self.config.get("record_section", "margins"))
            rectFactory = matplotlib_extras.axes.RectFactory(figure, margins=margins)
            ax = pyplot.axes(rectFactory.rect())
                
            vscale = self.config.getfloat("record_section", "%s_per_km" % field)

            for trace in streamC.traces:
                distKm = 1.0e-3 * (trace.stats.distance**2 + hypocenter.depth**2)**0.5
                t = trace.times(reftime=hypocenter.time)
                ax.plot(t, distKm*(1.0 + distKm*trace.data/vscale), linewidth=0.5, color="c_blue")
                if self.config.get("record_section", "show_labels"):
                    label = "%s.%s" % (trace.stats.network, trace.stats.station)
                    ax.text(numpy.min(t), distKm, label, horizontalalignment="left", verticalalignment="bottom", fontsize=6, color="c_orange")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Hypocentral Distance (km)")
            ax.set_xlim(0.0, tlength)
            ax.autoscale(axis="y", tight=True)

            if self.config.getboolean("record_section", "show_wavespeeds"):
                vp = self.config.getfloat("record_section", "vp_kmps")
                vs = self.config.getfloat("record_section", "vs_kmps")
                ymax = ax.get_ylim()[1]
                ax.plot([0.0, ymax/vp], [0.0, ymax], color="c_orange", linewidth=2, alpha=0.5)
                ax.text(ymax/vp, 0.995*ymax, "Vp={:3.1f}km/s".format(vp), color="c_orange", ha="left", va="top")
                ax.plot([0.0, ymax/vs], [0.0, ymax], color="c_red", linewidth=2, alpha=0.5)
                ax.text(ymax/vs, 0.995*ymax, "Vs={:3.1f}km/s".format(vs), color="c_red", ha="left", va="top")
            
            figure.savefig(os.path.join(plotsDir, "recordsection_%s_%s.pdf" % (field, component)))
        

        # Velocity
        return
    

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
    parser.add_argument("--plot-waveforms-map", action="store_true", dest="plot_waveforms_map")
    parser.add_argument("--plot-stations-map", action="store_true", dest="plot_stations_map")
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

    if args.plot_stations_map or args.all:
        plotter = StationsMap(app.params, args.show_progress)
        data.load_event()
        data.load_stations()
        plotter.plot(data)
    
    if args.plot_noise or args.all:
        plotter = NoiseFigure(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
        
    if args.plot_waveforms_map or args.all:
        plotter = WaveformsMap(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
    
    if args.plot_section or args.all:
        plotter = RecordSection(app.params, args.show_progress)
        data.load_event()
        data.load_processed_waveforms()
        plotter.plot(data)
    

# End of file
