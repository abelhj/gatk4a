FROM openjdk:8-jdk-slim
LABEL maintainer "Haley Abel <abelhj@wustl.edu>"

WORKDIR /opt/hall-lab

RUN apt-get update && \
        apt-get upgrade -y 


COPY build/libs/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar /opt/hall-lab/gatk-package-4.1.8.1-18-ge2f02f1-SNAPSHOT-local.jar