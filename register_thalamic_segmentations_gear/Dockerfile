# This makes a container for the gear. Replace the first path next to the second COPY function 
# below with a path to your local compiled matlab function

FROM neurodebian:xenial
MAINTAINER Ozzy(ozenctaskin@hotmail.com)

#Initialize flywheel v0 and move the required files into the container 
ENV FLYWHEEL /flywheel/v0/
RUN mkdir -p ${FLYWHEEL}
COPY manifest.json run ${FLYWHEEL}
RUN chmod +x /flywheel/v0/run

RUN apt-get install ants zip unzip python3 python3-pip
RUN pyton3-pip install nibabel==2.5.1 matplotlib imageio

# Set the entrypoint  
ENTRYPOINT /flywheel/v0/run
