#!/bin/bash
echo "--------------------------------------------------------------------"
echo "This script installs all necessary packages and compiles"
echo "Please refer to the script source code to see the list of packages"
echo "--------------------------------------------------------------------"

read -p "Please type 'y' to go ahead, any other key to exit: " -n 1 -r
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
	echo
	echo "Exiting."
	exit
fi

pacman -S --noconfirm --needed git
echo -e "\nInstalling compilation packages for building\n"
pacman -S --noconfirm --needed wget ${MINGW_PACKAGE_PREFIX}-cmake ${MINGW_PACKAGE_PREFIX}-gcc ${MINGW_PACKAGE_PREFIX}-make

./compile_on_mingw64.sh
