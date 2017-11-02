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

[processing.noise]
remove_bg = True
zero_coarse_levels = 0
preevent_threshold_reduction = 2.0
store_noise = False
store_orig = False
wavelet = coif4

[processing.filter]
zero_phase = True
num_corners = 2
freq_min = 0.2
freq_max = 25.0

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

[noise_figure]
height = 7.5
width = 20.0
margins = ((0.6, 0.5, 0.2), (0.5, 0.6, 0.7))


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
plots = plots
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
            if inventory is None:
                inventory = inventoryDC
            else:
                inventory += inventoryDC

        stream = self._get_raw_stream()
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
        storeNoise = self.params.getboolean("processing.noise", "store_noise")
        storeOrig = self.params.getboolean("processing.noise", "store_orig")
        accSM = obspyutils.noise.denoise(sSM, wavelet, removeBg, zeroCoarse, preWindow, preThreshold, storeOrig, storeNoise)

        freqMin = self.params.getfloat("processing.filter", "freq_min")
        freqMax = self.params.getfloat("processing.filter", "freq_max")
        numCorners = self.params.getint("processing.filter", "num_corners")
        zeroPhase = self.params.getboolean("processing.filter", "zero_phase")
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
        velBB = obspyutils.noise.denoise(sVel, wavelet, removeBg, zeroCoarse, preWindow, preThreshold, storeOrig, storeNoise)
        
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
        import obspyutils.baseline
        from basemap.Figure import Figure
        
        if self.showProgress:
            sys.stdout.write("Plotting noise figures...")
        
        self.figure = Figure()
        w = self.params.getfloat("noise_figure", "width")
        h = self.params.getfloat("noise_figure", "height")
        margins = literal_eval(self.params.get("noise_figure", "margins"))
        self.figure.open(w, h, margins)

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
        
        self.axes = {}
        for irow,row in enumerate(NoiseFigure.ROWS):
            for icol,col in enumerate(NoiseFigure.COLS):
                ax = self.figure.axes(nrows, ncols, irow+1, icol+1)
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
                self.axes["%s_%s" % (row,col)] = (ax,line)
                

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
        plotter = NoiseFigure(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
        
    if args.plot_map or args.all:
        plotter = WaveformsMap(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
    
    if args.plot_section or args.all:
        plotter = RecordSection(app.params, args.show_progress)
        data.load_processed_waveforms()
        plotter.plot(data)
    

# End of file
