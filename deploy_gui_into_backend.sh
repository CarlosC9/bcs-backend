#!/bin/bash
# Script to compile "bcs-gui" and deploy it into "bcs-backend", inside "static_gui" directory

cd ~/Proyectos/NEXTGENDEM/bcs-gui/
# Update Javascript packages
npm install
# Clear everything in the "dist" directory
rm dist -fr
# Compile
# ng build --prod --aot -c production_local --base-href /nis_client/
ng build --prod --aot --base-href /gui/

# Delete GUI files in the "bcs-backend" project
rm ~/Proyectos/NEXTGENDEM/bcs-backend/biobarcoding/static_gui/* -fr
# Copy just compiled bcs-gui files into "bcs-backend" project
cp -r ~/Proyectos/NEXTGENDEM/bcs-gui/dist/bcs-gui/* ~/Proyectos/NEXTGENDEM/bcs-backend/biobarcoding/static_gui
