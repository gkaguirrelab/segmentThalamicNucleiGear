# This makes a container for the gear. Replace the first path next to the second COPY function 
# below with a path to your local compiled matlab function

FROM freesurfer/freesurfer:7.1.1
MAINTAINER Ozzy(ozenctaskin@hotmail.com)

#Initialize flywheel v0 and move the required files into the container 
ENV FLYWHEEL /flywheel/v0/
RUN mkdir -p ${FLYWHEEL}
COPY manifest.json run ${FLYWHEEL}
RUN chmod +x /flywheel/v0/run

RUN (curl -sL https://rpm.nodesource.com/setup_12.x | bash -) \
  && yum clean all -y \
  && yum update -y \
  && yum install -y zip unzip nodejs tree libXt libXext ncurses-compat-libs python3 python3-pip \
  && yum clean all -y \
  && npm install npm --global

RUN source $FREESURFER_HOME/SetUpFreeSurfer.sh

# extra segmentations require matlab compiled runtime
RUN fs_install_mcr R2014b

# Delete the original tahalamic script to make sure it is gone
RUN rm $FREESURFER_HOME/bin/quantifyThalamicNuclei.sh

# Now replace it with the new thalamic script
COPY quantifyThalamicNuclei.sh $FREESURFER_HOME/bin/

# Set the entrypoint  
ENTRYPOINT /flywheel/v0/run
