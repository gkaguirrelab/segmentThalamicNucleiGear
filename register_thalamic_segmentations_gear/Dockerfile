# This makes a container for the gear. Replace the first path next to the second COPY function 
# below with a path to your local compiled matlab function

FROM neurodebian:bionic
MAINTAINER Ozzy(ozenctaskin@hotmail.com)

#Initialize flywheel v0 and move the required files into the container 
ENV FLYWHEEL /flywheel/v0/
RUN mkdir -p ${FLYWHEEL}
COPY manifest.json run ${FLYWHEEL}
RUN chmod +x /flywheel/v0/run
COPY fslinstaller.py /fslinstaller.py

# Install required packages    
RUN apt-get update \
    && apt-get install -y \
    git \
    python3 \
    python3-pip \
    zip \
    unzip \
    wget \ 
    mrtrix

RUN pip3 install --upgrade pip
RUN pip3 install 'nibabel==2.5.1' 'matplotlib' 'imageio'

# Install python2 for fslinstallation and install the newest verison of fsl
RUN apt-get install python -y
RUN python2 /fslinstaller.py -d /fsl 

# Download and install Freesurfer on call
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz
RUN mkdir freesurfer
RUN tar -zxf freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz 
COPY license.txt /freesurfer/license.txt
ENV FREESURFER_HOME /freesurfer

# Set the entrypoint  
ENTRYPOINT /flywheel/v0/run
