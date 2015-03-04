#!/usr/bin/env bash
#Install Script for corpus callosum thickness package for Linux/OSX
#
#__author__ = "Shantanu H. Joshi, Brandon R. Ayers"
#__copyright__ = "Copyright 2014 Shantanu H. Joshi 
#				  Ahmanson Lovelace Brain Mapping Center, University of California Los Angeles"
#__email__ = "sjoshi@bmap.ucla.edu"

message_if_failed() {
    #print out message and quit if last command failed
    if [ $? -ne 0 ]; then 
        echo -e >&2 $1
        exit 2
    fi
}

if [ $# -eq 0 ]
  then
    echo "usage: $0 <install directory>"
    echo "This program installs the corpus callosum thickness package in the specified install directory."
    exit 0
fi

orig_dir=`pwd`
install_dir=$1
mkdir $install_dir
mkdir $install_dir/"tmp"

osname=`uname`
if [[ "$osname" == "Darw"* ]]; then
	platform="MacOSX"
else
	platform="Linux"
fi;

R_data_table_package="data.table"
echo "'${R_data_table_package}' %in% rownames(installed.packages())" | R --quiet --no-save | grep -q "TRUE" #this will have exit code zero if package is installed
package_is_installed=$?
if [ $package_is_installed == 0 ]; then
    echo -e "\nGood, R package ${R_data_table_package} is already installed.\n"
else
	#Create a library path         
    R_lib_path=${install_dir}/Rlibs		
    mkdir $R_lib_path
    echo -e "\nR package ${R_data_table_package} is not installed, installing it now in ${R_lib_path}...\n"
    #install packages     
    echo "install.packages('${R_data_table_package}', repos='http://cran.us.r-project.org', lib='${R_lib_path}')" | R --quiet --no-save 1>> ${install_dir}/tmp/install.log
    message_if_failed "ERROR: failed to install data.table package for R. Exitting..."
    #Create a .Rprofile in user's HOME to make sure R can find this library.
    echo -e ".libPaths( c( .libPaths(), '${R_lib_path}'))" >> "${HOME}/.Rprofile" 
    message_if_failed "\nWARNING: Was not able to create a .Rprofile that would add ${R_lib_path} to the R's .libPaths()"
fi


R_lme4_package="lme4"
echo "'${R_lme4_package}' %in% rownames(installed.packages())" | R --quiet --no-save | grep -q "TRUE" #this will have exit code zero if package is installed
package_is_installed=$?
if [ $package_is_installed == 0 ]; then
    echo -e "\nGood, R package ${R_lme4_package} is already installed.\n"
else
	#Create a library path         
    R_lib_path=${install_dir}/Rlibs		
    mkdir $R_lib_path
    echo -e "\nR package ${R_lme4_package} is not installed, installing it now in ${R_lib_path}...\n"
    #install packages     
    echo "install.packages('${R_lme4_package}', repos='http://cran.us.r-project.org', lib='${R_lib_path}')" | R --quiet --no-save 1>> ${install_dir}/tmp/install.log
    message_if_failed "ERROR: failed to install data.table package for R. Exitting..."
    #Create a .Rprofile in user's HOME to make sure R can find this library.
    echo -e ".libPaths( c( .libPaths(), '${R_lib_path}'))" >> "${HOME}/.Rprofile" 
    message_if_failed "\nWARNING: Was not able to create a .Rprofile that would add ${R_lib_path} to the R's .libPaths()"
fi

echo "Downloading anaconda python...This may take a few minutes..."
curl -o $install_dir/tmp/Miniconda-3.4.2-${platform}-x86_64.sh http://repo.continuum.io/miniconda/Miniconda-3.4.2-${platform}-x86_64.sh
chmod +x $install_dir/tmp/Miniconda-3.4.2-${platform}-x86_64.sh
echo "Done."
echo -n "Installing anaconda python..."
$install_dir/tmp/Miniconda-3.4.2-${platform}-x86_64.sh -b -f -p $install_dir 1> ${install_dir}/tmp/install.log
echo "Done."
echo -n "Installing statsmodels...This may take a few minutes..."
$install_dir/bin/conda install statsmodels -q --yes 1>> ${install_dir}/tmp/install.log
$install_dir/bin/conda install pip -q --yes 1>> ${install_dir}/tmp/install.log
echo "Done."

echo -n "Installing vtk..."
$install_dir/bin/conda install vtk -q --yes 1>> ${install_dir}/tmp/install.log
echo -n "Installing rpy2..."
#$install_dir/bin/pip install -q rpy2 1>> ${install_dir}/tmp/install.log
$install_dir/bin/pip install rpy2==2.3.9 1>> ${install_dir}/tmp/install.log
echo "Done."

echo "R_LIBS=${R_lib_path}:$R_LIBS; export R_LIBS" >> ~/.bashrc
echo "R_LIBS=${R_lib_path}:$R_LIBS; export R_LIBS" >> ~/.bash_profile
echo "setenv R_LIBS ${R_lib_path}:$R_LIBS" >> ~/.tcshrc
source ~/.bashrc

echo -n "Installing shapestats package...This may take a few minutes..."
$install_dir/bin/pip install git+https://github.com/bmapdev/shapestats
echo -n "Installing shapeio package...This may take a few minutes..."
$install_dir/bin/pip install git+https://github.com/bmapdev/shapeio
echo -n "Installing curvematch package...This may take a few minutes..."
$install_dir/bin/pip install git+https://github.com/bmapdev/curvematch
echo -n "Installing corpus callosum shape package...This may take a few minutes..."
$install_dir/bin/pip install git+https://github.com/bmapdev/ccshape
echo "Done."

echo "Install complete. Clearing package storage..."
rm -r $install_dir/pkgs/
rm -r $install_dir/tmp/
echo "Corpus callosum thickness package was installed successfully"
exit 0