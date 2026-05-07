FROM ubuntu:noble

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        python3 \
        python3-pip \
        python3-dev \
        bc \
        time \
        ca-certificates \
        git \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --break-system-packages \
        biopython \
        bcbio-gff \
        matplotlib

WORKDIR /opt/panprova

COPY . /opt/panprova

RUN bash compile.sh

ENV PANPROVA_PATH=/opt/panprova
ENV PATH=/opt/panprova:${PATH}

WORKDIR /work

CMD ["bash"]
