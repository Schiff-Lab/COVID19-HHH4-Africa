#!/bin/sh
cd path_to_folder_with_R_scripts
Rscript 01.R
Rscritp 02.R
Rscript 03.R
git add -A
git commit -m "Updated model with new data"
git push
sudo sytemctl stop shiny-server
systemctl start shiny-server
echo mypassword | sudo -S command
echo "Task has been run :)"
