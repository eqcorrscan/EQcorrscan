"""
Script to validate the synthetic detector method against test earthquakes
"""

def synth_from_sfile(sfile):
    """
    Function to generate a synthetic template for a given s-file

    :type path: str
    :param path: Path to the sfile
    """
    from eqcorrscan.utils import Sfile_util
    from eqcorrscan.utils import syth_seis
    picks=Sfile_util.readpicks(sfile)
    # We only want P and S phases
    picks=[p for p in picks if p.phase in ['P','S']]
    # We want a list of the stations that we have picks for
    stations=list(set([p.station for p in picks]))
    # Loop through the stations
    synths=[]
    for station in stations:
        
