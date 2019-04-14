FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
    apt-get install -y wget curl unzip python3 python3-pip bats \
    awscli libcurl4-openssl-dev libhdf5-dev python-tables hdf5-tools

RUN pip3 install pandas>=0.22.0 boto3>=1.7.2 feather-format \
                 s3fs tables scipy joblib==0.13.2

# Add the script to the PATH
ADD ./make-experiment-collection.py /usr/local/bin/
ADD lib /usr/local/bin/lib

RUN mkdir /scratch

# Run tests
ADD tests/ /usr/local/tests
RUN bats /usr/local/tests && \
    rm -r /usr/local/tests
