#!/bin/bash
#Add the direct commands for these scripts to your .bashrc file.

here=$(pwd)

echo ' ' >> ~/.bashrc
echo \#\#\# Spectroscopy commands: >> ~/.bashrc
echo alias degrade=$here/'degrade.py' >> ~/.bashrc
echo alias seehead=$here/'seehead.py' >> ~/.bashrc
echo alias seefits=$here/'plotSpec.py' >> ~/.bashrc
echo alias cutSpec=$here/'cutSpec.py' >> ~/.bashrc
echo alias fits2txt=$here/'fits2txt.py' >> ~/.bashrc
echo alias txt2fits=$here/'txt2fits.py' >> ~/.bashrc
echo \#\#\# >> ~/.bashrc
echo ' ' >> ~/.bashrc
