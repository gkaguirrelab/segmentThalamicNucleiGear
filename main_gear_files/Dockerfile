# This makes a container for the gear. Replace the first path next to the second COPY function 
# below with a path to your local compiled matlab function

FROM freesurfer/freesurfer:7.1.1
MAINTAINER Ozzy(ozenctaskin@hotmail.com)

#Initialize flywheel v0 and move the required files into the container 
ENV FLYWHEEL /flywheel/v0/
RUN mkdir -p ${FLYWHEEL}
COPY manifest.json run ${FLYWHEEL}
COPY segmentTocifti_files.zip /tmp/segmentTocifti_files.zip
COPY xfm_subcortToHcp.py ${FLYWHEEL}/xfm_subcortToHcp.py
RUN chmod +x /flywheel/v0/run
RUN chmod +x ${FLYWHEEL}/xfm_subcortToHcp.py

RUN (curl -sL https://rpm.nodesource.com/setup_12.x | bash -) \
  && yum clean all -y \
  && yum update -y \
  && yum install -y zip unzip nodejs tree libXt libXext ncurses-compat-libs python3 python3-pip git openssl git \
  && yum groupinstall -y "Development Tools" \
  && yum clean all -y \
  && npm install npm --global

RUN source $FREESURFER_HOME/SetUpFreeSurfer.sh

# extra segmentations require matlab compiled runtime
RUN fs_install_mcr R2014b

# Delete the original tahalamic script to make sure it is gone
RUN rm $FREESURFER_HOME/bin/quantifyThalamicNuclei.sh

# Now replace it with the new thalamic script
COPY quantifyThalamicNuclei.sh $FREESURFER_HOME/bin/

# Install workbench 
RUN yum search freetype && yum install freetype -y
RUN wget https://www.humanconnectome.org/storage/app/media/workbench/workbench-rh_linux64-v1.5.0.zip
RUN unzip workbench-rh_linux64-v1.5.0.zip
RUN cp workbench/bin_rh_linux64/* /usr/bin; cp -r workbench/* /usr/

# Install ants
RUN yum install -y openssl-devel
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2.tar.gz
RUN tar -xvf cmake-3.20.2.tar.gz
RUN cd cmake-3.20.2; ./configure; gmake
RUN cp -r cmake-3.20.2/bin/* /usr/bin
RUN cp -r cmake-3.20.2 /usr/share/cmake-3.20
RUN git clone https://github.com/cookpa/antsInstallExample.git
RUN cd antsInstallExample; bash installANTs.sh
RUN cp antsInstallExample/install/bin/* /usr/bin/
RUN cp antsInstallExample/install/lib/* /usr/lib/

# Set the entrypoint  
ENTRYPOINT /flywheel/v0/run
