---
published: false
---
## How to Install RStudio Server on WSL - Ubuntu 18.04

If you want to have RStudio Server on your computer running Windows 10, for whatever reason, this is a short tutorial on how to do that. One of my motivations is that brms models not always run successfully on Windows.

I assume you already have installed WSL and have some grasp on how Linux works.

First you need to add a GPG key and the R repository to apt.

```
curl -sL "http://keyserver.ubuntu.com/pks/lookup?op=get&search=0x51716619E084DAB9" | sudo apt-key add
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
```

Then you have to update your local apt library, the packages and install R and gdebi (used for installing software packages in .deb format)
```
sudo apt update
sudo apt upgrade
sudo apt install r-base gdebi-core
```

Then you have to download the correct version of libssl (a dependency for RStudio Server) and install it.
```
wget http://ftp.us.debian.org/debian/pool/main/o/openssl1.0/libssl1.0.2_1.0.2q-2_amd64.deb
sudo gdebi libssl1.0.2_1.0.2q-2_amd64.deb
```

After that step, you could install RStudio Server, and then start it. You have to start it after every reboot.
```
wget https://download2.rstudio.org/rstudio-server-1.1.463-amd64.deb
sudo gdebi rstudio-server-1.1.463-amd64.deb
sudo rstudio-server start
```

Then you go to localhost:8787 to access RStudio Server.

Voila!