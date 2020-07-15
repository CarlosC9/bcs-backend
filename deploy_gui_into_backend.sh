#!/bin/bash
# Script to compile "bcs-gui" and deploy it into "bcs-backend", inside "static_gui" directory

cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-gui/
# Update Javascript packages
npm install
# Clear everything in the "dist" directory
rm dist -fr
# Compile
# ng build --prod --aot -c production_local --base-href /nis_client/
ng build --prod --aot --base-href /gui/

# Delete GUI files in the "bcs-backend" project
rm ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/biobarcoding/static_gui/* -fr
# Copy just compiled bcs-gui files into "bcs-backend" project
cp -r ~/GoogleDrive/AA_NEXTGENDEM/bcs-gui/dist/bcs-gui/* ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/biobarcoding/static_gui
