## Container Setup

To run the dashboard, ensure you have moved into the directory where the Dockerfile is located, and run the command below in a terminal:
```
docker build -t geotiled-dashboard .
```
```
docker run --rm -dp 10142:10142 geotiled-dashboard
```
From there, you should be able to access the dashboard at [http://localhost:10142/dashboard](http://localhost:10142/dashboard).
    > If running from a VM, instead replace `localhost` with the address of the VM.