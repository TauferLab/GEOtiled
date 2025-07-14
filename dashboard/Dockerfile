# python slim image 
FROM python:3.10-slim

RUN apt update && apt install -y \ 
    libgdal-dev=3.6.2+dfsg-1+b2

RUN apt-get update && apt-get install -y \
    wget \
    g++ \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN pip install boto3==1.35.99 python-dotenv pyyaml numpy geopandas matplotlib requests pandas tqdm ipywidgets ipywidgets_bokeh ipykernel watchfiles panel==1.6.0 wheel setuptools>=67

RUN pip install gdal[numpy]==3.6.2

WORKDIR /usr/src/dashboard
COPY ./dashboard.py ./dashboard.py
COPY ./tools.py ./tools.py
COPY ./config.yaml ./config.yaml
COPY ./.env ./.env
EXPOSE 10142

ENV BOKEH_ALLOW_WS_ORIGIN="*"

CMD ["python", "-m", "panel", "serve", "dashboard.py", "--address", "0.0.0.0", "--allow-websocket-origin", "*", "--port", "10142"]
