# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure("2") do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://vagrantcloud.com/search.
  config.vm.box = "sylabs/singularity-3.6-ubuntu-bionic64"

  # Increase the default disk size for this Vagrant box from 20 GB to 30 GB:
  # config.vm.disk :disk, size: "30GB", name: "extra_storage", primary: true
  config.disksize.size = '20GB'

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  # NOTE: This will enable public access to the opened port
  # config.vm.network "forwarded_port", guest: 80, host: 8080

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine and only allow access
  # via 127.0.0.1 to disable public access
  # config.vm.network "forwarded_port", guest: 80, host: 8080, host_ip: "127.0.0.1"

  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
  # config.vm.network "private_network", ip: "192.168.33.10"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  # config.vm.synced_folder "../data", "/vagrant_data"
  #config.vm.synced_folder "/Users/chrisjackson/Desktop/VMs/vm-singularity/data", "/vagrant_data"
  

  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  config.vm.provider "virtualbox" do |vb|
  #   # Display the VirtualBox GUI when booting the machine
  #   vb.gui = true
  #
    # Customize the amount of memory on the VM:
    vb.memory = "5120"
    vb.cpus = "4"
    vb.name = "hybpiper-paragone-vm"
  end
  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Enable provisioning with a shell script. Additional provisioners such as
  # Ansible, Chef, Docker, Puppet and Salt are also available. Please see the
  # documentation for more information about their specific syntax and use.
  config.vm.provision "shell", inline: <<-SHELL

    # Re-size disk:
    parted /dev/sda resizepart 1 100%
    pvresize /dev/sda1
    lvresize -rl +100%FREE /dev/mapper/vagrant--vg-root

    
    # Install general packages:
    apt-get install -y software-properties-common
    apt-get update
    apt-get install -y curl
    apt install -y vim
    apt install -y screen
    apt install -y default-jre
    apt install -y unzip


    # Install R version 4:
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
    apt update
    apt install -y r-base r-base-core r-recommended r-base-dev


    # Install R packages:
    add-apt-repository ppa:c2d4u.team/c2d4u4.0+
    apt update
    apt install -y --no-install-recommends r-cran-ape r-cran-stringr r-cran-seqinr


    # Create bin directory in /home/vagrant:
    runuser -l vagrant -c 'mkdir bin' 


    # Install miniconda in /home/vagrant/bin:
    runuser -l vagrant -c '
    if [ ! -d /home/vagrant/bin/miniconda3 ]; then \
      curl -OL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
      bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/vagrant/bin/miniconda3 && \
      rm Miniconda3-latest-Linux-x86_64.sh; \
    fi'
  

    # Add conda to PATH:
    echo 'export PATH="$PATH:/home/vagrant/bin/miniconda3/bin"' >> /home/vagrant/.bashrc


    # Add bioconda channel:
    runuser -l vagrant -c '
    /home/vagrant/bin/miniconda3/bin/conda config --add channels defaults; \
    /home/vagrant/bin/miniconda3/bin/conda config --add channels bioconda; \
    /home/vagrant/bin/miniconda3/bin/conda config --add channels conda-forge'


    # Install programs using conda:
    runuser -l vagrant -c '
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::iqtree=2.1.2; \
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::mafft=7.475; \
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::bwa=0.7.17; \
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::samtools=1.9; \
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::bbmap=38.86; \
    /home/vagrant/bin/miniconda3/bin/conda install -y bioconda::bcftools=1.9'


    # Install Astral in /home/vagrant/bin:
    runuser -l vagrant -c '
    cd bin; \
    curl -OL https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.7.zip; \
    unzip Astral.5.7.7.zip; \
    rm Astral.5.7.7.zip; \
    chmod -R a+rX Astral; \
    cd ..'


    # Install Nextflow: 
    runuser -l vagrant -c '
    cd bin; \
    mkdir nextflow_install; \
    cd nextflow_install; \ 
    curl -s https://get.nextflow.io | bash; \
    cd ../..'


    # Add Nextflow to PATH:
    echo 'export PATH="$PATH:/home/vagrant/bin/nextflow_install"' >> /home/vagrant/.bashrc


    # Clone the HybPiper-RBGV repo:
    runuser -l vagrant -c '
    cd bin; \
    git clone https://github.com/chrisjackson-pellicle/hybpiper-nf.git; \
    cd ..'


    # Clone the Yang-and-Smith repo:
    runuser -l vagrant -c '
    cd bin; \
    git clone https://github.com/chrisjackson-pellicle/paragone-nf.git; \
    cd ..'


    # Add screen tab captioning:
    echo 'caption always "%{= kw}%-w%{= BW}%n %t%{-}%+w %-= @%H - %LD %d %LM - %c"' > /home/vagrant/.screenrc

  SHELL
end
