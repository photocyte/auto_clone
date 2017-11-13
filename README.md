# auto_clone

This is a standalone python based server that intends to make the production of gibson primers & in-silico plasmids more straightforward.

Dependencies:
* Python3.6+
* BioPython
* Pandas
* SciPy (for the Tornado stand alone webserver)
* xlwt (for export of Excel file)

Installation (using [MacPorts](https://www.macports.org)):
```
git clone https://github.com/photocyte/auto_clone.git
sudo port install py36
sudo port install py36-scipy
sudo port install py36-pandas
sudo port install py36-biopython
sudo port install py36-xlwt
```

Usage:
```
sudo port select --set python3 python36
cd auto_clone
python3 ./auto_clone_server.py
```
Default webaddress is http://127.0.0.1:5686

See [how_to_add_plasmids.txt](./how_to_add_plasmids.txt) to add plasmids.


