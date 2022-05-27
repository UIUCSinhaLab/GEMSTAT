FROM ubuntu:20.04 as gemstat_build
LABEL maintainer="bryan.j.lunt+gemstat@gmail.com"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get install -y curl && \
    apt-get install -y autoconf automake libtool pkg-config && \
    apt-get install -y libtool libtool-bin libltdl7 && \
    apt-get install -y gcc g++ && \
    apt-get install -y cmake && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y -q && \
    apt-get install -y git vim && \
    rm -rf /var/lib/apt/lists/*
    
RUN apt-get update -y -q && \
    apt-get install -y libgsl-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
