# auto_clone

This is a standalone python based server that intends to make the production of gibson primers & in-silico plasmids more straightforward.

Dependencies:
* Python3.6+
* BioPython
* Pandas
* SciPy (for the Tornado stand alone webserver)

Installation:
* sudo port install py36
* sudo port install py36-scipy
* sudo port install py36-pandas
* sudo port install py36-biopython

Usage:
* python auto_clone_server.py
* Default webaddress is http://127.0.0.1:5686
