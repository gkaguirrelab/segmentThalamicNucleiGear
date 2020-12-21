# This makes a container for the gear. Replace the first path next to the second COPY function 
# below with a path to your local compiled matlab function

FROM neurodebian:xenial
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
RUN pip3 install nibabel==2.5.1 matplotlib imageio

# Install python2 for fslinstallation and install the newest verison of fsl
RUN apt-get install python -y
RUN python2 /fslinstaller.py -d /fsl 

# Set the entrypoint  
ENTRYPOINT /flywheel/v0/run
