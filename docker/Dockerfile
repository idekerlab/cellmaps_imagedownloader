# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get --allow-releaseinfo-change update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/cellmaps_imagedownloader
COPY ./ /tmp/cellmaps_imagedownloader/
RUN pip install /tmp/cellmaps_imagedownloader

RUN rm -rf /tmp/cellmaps_imagedownloader

ENTRYPOINT ["/opt/conda/bin/cellmaps_imagedownloadercmd.py"]

CMD ["--help"]
