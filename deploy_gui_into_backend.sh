#!/bin/bash
# Script to compile "bcs-gui" and deploy it into "bcs-backend", inside "static_gui" directory

if ["$(whoami)" == "rnebot"]; then
  ngddir = ~/GoogleDrive/AA_NEXTGENDEM
elif [ "$(whoami)" == "acurbelo" ]; then
  ngddir = ~/Proyectos/NEXTGENDEM
elif [ "$(whoami)" == "pmoreno" ]; then
  ngddir = ~/pmoreno/NEXTGENDEM
elif [ "$(whoami)" == "dreyes" ]; then
  ngddir = ~/dreyes/NEXTGENDEM
fi

# CD
cd $ngddir/bcs-gui
# Update Javascript packages
npm install
# Clear everything in the "dist" directory
rm dist -fr
# Compile
# ng build --prod --aot -c production_local --base-href /nis_client/
ng build --prod --aot --base-href /gui/

# Delete GUI files in the "bcs-backend" project
rm $ngddir/bcs-backend/biobarcoding/static_gui/* -fr
# Copy just compiled bcs-gui files into "bcs-backend" project
cp -r $ngddir/bcs-gui/dist/bcs-gui/* $ngddir/bcs-backend/biobarcoding/static_gui
