---
title: "Deploying BIGapp with Docker/singularity"
description:"Access BIGapp interface using BIGapp's Docker image stored at DockerHub"
date: 2025-08-22
categories: [installation]
tags: [bigapp, install]
layout: post
author: cris
mermaid: true
---

## What are Docker containers?

Docker containers bundle an application together with all its dependencies (system libraries, R packages, config) into a portable image. That means:

* Reproducibility: everyone runs the same environment—no "works on my machine."
* Fast setup: skip lengthy installs; just pull and run.
* Portability: run the same image on laptops, servers, or HPC via Apptainer/Singularity.
* Isolation: avoid conflicts with other software on the host.

## Overview

We maintain a Docker image for BIGapp that updates automatically whenever the package is updated on GitHub. Using the image ensures the app and all dependencies are exactly as intended, saving installation time and improving reproducibility. It also makes it easy to run BIGapp on HPC systems to speed up analyses.

Below are three ways to access BIGapp via containers.

## Local (Docker Desktop GUI)

1. Open Docker Desktop → Docker Hub (left sidebar).
2. In the search box, type: breedinginsight/bigapp and open it.
3. Go to the tag tab and one you want (e.g., latest or a version like v1.5.1).
4. Click Pull and wait for it to complete.

Now it will appear under Images:

1. Go to the Images tab and find breedinginsight/bigapp (the tag you pulled).
2. Click Run. In the “Run a new container” panel:
* Name: bigapp (or anything).
* Ports → Publish a new network port
  * Host port: 8080
  * Leave environment variables empty (the image serves on port 80).
  * Click Run.

3. Open it:

* In Containers, click your container, then click the 8080:80 port chip (or Open in Browser), or manually open http://localhost:8080.

Tip: If port 8080 is already in use, pick another host port (e.g., 9000) and map 9000:80.

## Local (Command Line: Docker/Podman)

```bash
# 1) Pull the image (pick a tag)
docker pull docker.io/breedinginsight/bigapp:latest

# 2) Run it (host port 8080 -> container port 80)
docker run -d --name bigapp \
  --restart unless-stopped \
    -p 8080:80 \
      docker.io/breedinginsight/bigapp:latest

# 3) Open in your browser:
# http://localhost:8080

# 4) Stop and remove the container when done
docker rm -f bigapp
```

If 8080 is taken: change the left side of -p, e.g. -p 9000:80 then visit http://localhost:9000.
Podman users: replace docker with podman (same flags).

## HPC (Apptainer/Singularity + Slurm, no sudo)

```bash

# 1) On the login node: pull the image once (creates a .sif file)
apptainer pull bigapp_latest.sif docker://breedinginsight/bigapp:latest

# 2) Request an interactive compute allocation (adjust for your cluster)
salloc -c 4 --mem=16G --time=02:00:00

# 3) On the compute node shell that opens, run the app on port 8080
export PORT=8080
apptainer run --cleanenv --env PORT=$PORT bigapp_latest.sif
```

The app now listens on localhost:8080 on the compute node.

From your laptop (new terminal), create an SSH tunnel to reach the compute node via the login node:

```bash
# Replace user and host with your cluster login node.
# If your cluster requires a direct tunnel to the compute node, adapt accordingly.
ssh -L 8080:localhost:8080 [userID]@cbsubi2.biohpc.cornell.edu
```

Now open http://localhost:8080 in your browser.

When finished:

* Press Ctrl+C in the SSH tunnel terminal to close the tunnel.
* Type exit in the compute node shell to end the Slurm job.

Notes
* Some clusters expose compute nodes differently. If your job runs on a node like `cn123`, you may need: `ssh -L 8080:cn123:8080 your_user@login.node` (or use `ssh -J` for a jump host).

* If your site uses `singularity` instead of `apptainer`, replace `apptainer` with `singularity`.
* Adjust `-c`, `--mem`, and `--time` to fit your workload and queue policies.
