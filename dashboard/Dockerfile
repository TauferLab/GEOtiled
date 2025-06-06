# build the conda environment and file structure
FROM continuumio/miniconda3:latest AS build

WORKDIR /usr/src/dashboard
COPY ./environment.yaml ./environment.yaml

RUN apt-get update && apt-get install -y \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN conda env create -f ./environment.yaml && \
    conda clean -a -y && \
    rm -rf /opt/conda/pkgs/* /root/.cache/* /var/lib/apt/lists/*

# python slim image 
FROM python:3.10-slim

ENV PATH=/opt/conda/bin:$PATH

COPY --from=build /opt/conda /opt/conda
COPY --from=build /usr/src/dashboard /usr/src/dashboard

WORKDIR /usr/src/dashboard
COPY ./dashboard.py ./dashboard.py
COPY ./tools.py ./tools.py
COPY ./config.yaml ./config.yaml
COPY ./.env ./.env
EXPOSE 10142

ENV BOKEH_ALLOW_WS_ORIGIN="*"

CMD ["conda", "run", "--no-capture-output", "-n", "geotiled", "python", "-m", "panel", "serve", "dashboard.py", "--address", "0.0.0.0", "--allow-websocket-origin", "*", "--port", "10142"]
