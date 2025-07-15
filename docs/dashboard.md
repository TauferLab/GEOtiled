# 📊 Dashboard Guide

This guide goes over how to create and run the dashboard for data are generated with GEOtiled. Currently, 17 different terrain parameters generated with GEOtiled at 30m and 10m resolution within the region of the Continental United States (CONUS) are available for visualization.

## 🚀 Quick Start

Before installation, ensure you have [Docker](https://www.docker.com/get-started/) installed.

### Build Container

After cloning the GEOtiled repository, move into the `dashboard` directory. Run the following in the terminal to build the container.

```
docker build -t geotiled-dashboard .
```

### Run container

Once the container is built, run the following in the terminal to activate the container.

```
docker run --rm -dp 10142:10142 geotiled-dashboard
```

From there, you should be able to access the dashboard at [http://localhost:10142/dashboard](http://localhost:10142/dashboard).
> If running from a VM, instead replace `localhost` with the address of the VM.

### Notes

* The larger the data, the longer it will take to visualize.
  * Larger states and higher resolutions directly correlate to larger file size.
* It is recommended to downsample larger data for performance improvements to visualization.